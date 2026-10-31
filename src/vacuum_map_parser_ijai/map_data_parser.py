"""Ijai map parser."""

import logging
import math
import zlib
from dataclasses import replace
from pathlib import Path as FilePath
from typing import Any

from PIL import ImageFont

from vacuum_map_parser_base.config.color import ColorsPalette, SupportedColor
from vacuum_map_parser_base.config.drawable import Drawable
from vacuum_map_parser_base.config.image_config import ImageConfig
from vacuum_map_parser_base.config.size import Size, Sizes
from vacuum_map_parser_base.config.text import Text
from vacuum_map_parser_base.image_generator import ImageGenerator
from vacuum_map_parser_base.map_data import Area, ImageData, MapData, Path, Point, Room, Wall, Zone
from vacuum_map_parser_base.map_data_parser import MapDataParser

import vacuum_map_parser_ijai.beautify_min as Beautify
import vacuum_map_parser_ijai.RobotMap_pb2 as RobotMap

from .aes_decryptor import decrypt
from .ijai_coordinate_transforms import Transformer
from .image_parser import IjaiImageParser

_LOGGER = logging.getLogger(__name__)

# Fraction of the final rendered image width the vacuum icon's radius
# should take up, matched by measuring the official Xiaomi app's icon
# against its map image.
TARGET_VACUUM_RADIUS_RATIO = 0.022

# Fraction of the final rendered image width used as the room-name font
# size.
ROOM_NAME_FONT_RATIO = 0.035
ROOM_NAME_MIN_FONT_SIZE = 10

# Bundled font used for room names: Pillow's built-in default font doesn't
# cover Turkish characters (ı, ğ, ş, ö, ü, ç render as missing-glyph boxes),
# and relying on whatever happens to be installed on the host isn't
# reliable -- Home Assistant OS ships with no system fonts.
_FONT_PATH = FilePath(__file__).parent / "fonts" / "DejaVuSans.ttf"


class _IjaiImageGenerator(ImageGenerator):
    """ImageGenerator that sizes room-name text relative to the map instead
    of using vacuum_map_parser_base's hardcoded-tiny default PIL font."""

    def _draw_room_names(self, map_data: MapData) -> None:
        if map_data.rooms is None or map_data.image is None:
            return
        color = self._get_color(SupportedColor.ROOM_NAMES)
        dims = map_data.image.dimensions
        font_size = max(ROOM_NAME_MIN_FONT_SIZE, int(dims.width * dims.scale * ROOM_NAME_FONT_RATIO))
        font: Any
        try:
            font = ImageFont.truetype(str(_FONT_PATH), font_size)
        except OSError:
            _LOGGER.warning("Unable to load bundled font %s, falling back to PIL default", _FONT_PATH)
            font = ImageFont.load_default(size=font_size)
        for room in map_data.rooms.values():
            p = room.point()
            if p is None or room.name is None:
                continue
            point = p.to_img(dims)
            self._draw_room_name_text(map_data.image, room.name, point.x, point.y, font, color)

    @staticmethod
    def _draw_room_name_text(image: ImageData, text: str, x: float, y: float, font: Any, color: Any) -> None:
        def draw_func(draw: Any) -> None:
            left, top, right, bottom = draw.textbbox((0, 0), text, font)
            w, h = right - left, bottom - top
            draw.text((x - w / 2, y - h / 2), text, font=font, fill=color)

        ImageGenerator._draw_on_new_layer(image, draw_func, ImageGenerator._use_transparency(color))


class IjaiMapDataParser(MapDataParser):
    """Ijai map parser."""

    POSITION_UNKNOWN = 1100
    VIRTUALWALL_TYPE_WALL = 2
    VIRTUALWALL_TYPE_NO_MOP = 6
    VIRTUALWALL_TYPE_NO_GO = 3

    # pylint: disable=E1101
    robot_map = RobotMap.RobotMap()

    def __init__(
        self,
        palette: ColorsPalette,
        sizes: Sizes,
        drawables: list[Drawable],
        image_config: ImageConfig,
        texts: list[Text]
    ):
        super().__init__(palette, sizes, drawables, image_config, texts)
        self._image_parser = IjaiImageParser(palette, image_config, drawables)
        self._drawables = drawables
        self._image_generator = _IjaiImageGenerator(palette, sizes, drawables, image_config, texts)

    def unpack_map(self, raw_encoded: bytes, *args: Any, **kwargs: Any) -> bytes:
        return zlib.decompress(
            decrypt(
                raw_encoded,
                kwargs['wifi_sn'],
                kwargs['owner_id'],
                kwargs['device_id'],
                kwargs['model'],
                kwargs['device_mac']))

    def parse(self, raw: bytes, *args: Any, **kwargs: Any) -> MapData:
        map_data = MapData(0, 1)

        self.robot_map.ParseFromString(raw)
        # pylint: disable=W0201
        self.coord_transformer = Transformer(self.robot_map)

        if hasattr(self.robot_map, "mapData"):
            map_data.image, map_data.rooms, map_data.cleaned_rooms = self._parse_image()

        if hasattr(self.robot_map, "historyPose"):
            map_data.path = IjaiMapDataParser._parse_history()

        if hasattr(self.robot_map, "chargeStation"):
            pos_info = self.robot_map.chargeStation
            map_data.charger = Point(
                x=pos_info.x, y=pos_info.y, a=pos_info.phi * 180 / math.pi)
            _LOGGER.debug("pos: %s", map_data.charger)

        if hasattr(self.robot_map, "currentPose"):
            pos_info = self.robot_map.currentPose
            map_data.vacuum_position = Point(
                x=pos_info.x, y=pos_info.y, a=pos_info.phi * 180 / math.pi)
            _LOGGER.debug("pos: %s", map_data.vacuum_position)

        if (
                hasattr(self.robot_map, "mapInfo")
                and hasattr(self.robot_map, "roomDataInfo")
                and map_data.rooms is not None):
            IjaiMapDataParser._parse_rooms(map_data.rooms)

        if hasattr(self.robot_map, "virtualWalls"):
            (map_data.walls,
             map_data.no_go_areas,
             map_data.no_mopping_areas) = IjaiMapDataParser._parse_restricted_areas()

        if hasattr(self.robot_map, "areasInfo"):
            map_data.zones = IjaiMapDataParser._parse_cleaning_zones()

        if hasattr(self.robot_map, "navigationPoints"):
            map_data.goto = IjaiMapDataParser._parse_goto_point()

        if map_data.rooms is not None:
            _LOGGER.debug("rooms: %s", [str(room)
                          for number, room in map_data.rooms.items()])
            if map_data.rooms is not None and len(map_data.rooms) > 0 and map_data.vacuum_position is not None:
                vacuum_position_on_image = self.coord_transformer.map_to_image(
                    map_data.vacuum_position)
                map_data.vacuum_room = IjaiImageParser.get_current_vacuum_room(
                    self.robot_map.mapData.mapData, vacuum_position_on_image, IjaiMapDataParser.robot_map.mapHead.sizeX)
                if map_data.vacuum_room is not None:
                    map_data.vacuum_room_name = map_data.rooms[map_data.vacuum_room].name
                _LOGGER.debug("current vacuum room: %s", map_data.vacuum_room)

        if map_data.image is not None and not map_data.image.is_empty:
            # Icon/line sizes (vacuum, charger, obstacles, path width, ...)
            # are plain pixel values that vacuum_map_parser_base never
            # scales, so they end up a fixed, inconsistent fraction of the
            # final image depending on how much the auto-crop/auto-scale in
            # IjaiImageParser enlarged this particular map. Instead of
            # scaling them relative to the map's own (highly variable)
            # scale factor, target a fixed fraction of the final rendered
            # width -- matching how big the vacuum icon looks in the
            # official app regardless of how tightly this map got cropped.
            dims = map_data.image.dimensions
            final_width = dims.width * dims.scale
            base_vacuum_radius = self._sizes.get_size(Size.VACUUM_RADIUS)
            target_vacuum_radius = TARGET_VACUUM_RADIUS_RATIO * final_width
            icon_scale = target_vacuum_radius / base_vacuum_radius if base_vacuum_radius else 1
            if icon_scale != 1:
                scaled_sizes = Sizes({size: self._sizes.get_size(size) * icon_scale for size in Size})
                self._image_generator = _IjaiImageGenerator(
                    self._palette, scaled_sizes, self._drawables, self._image_config, self._texts)
            self._image_generator.draw_map(map_data)

        return map_data

    def _parse_image(self) -> tuple[ImageData, dict[int, Room], set[int]]:
        image_left = 0
        image_top = 0
        image_width = self.robot_map.mapHead.sizeX
        image_height = self.robot_map.mapHead.sizeY
        image_size = image_height * image_width
        _LOGGER.debug("width: %d, height: %d", image_width, image_height)

        # Non painted map tranformation
        will_beautify = (
            len(set(self.robot_map.mapData.mapData).symmetric_difference(
                [0, 128, 127])) == 0
            and len(self.robot_map.roomChain) > 0
            and self.robot_map.mapType == 0)
        if will_beautify:
            buautify_obj = Beautify.BeautifyMap(self.robot_map.mapHead)
            buautify_obj.setMap(self.robot_map.mapData)
            buautify_obj.transform()
            buautify_obj.roomColorByChain(self.robot_map.roomChain)
            buautify_obj.fillInternalObstacles()
            buautify_obj.normalizeMap()
            self.robot_map.mapData.mapData = bytes(buautify_obj.getMap())

        image, rooms_raw, cleaned_areas, effective_trim, effective_scale = self._image_parser.parse(
            self.robot_map.mapData.mapData, image_width, image_height)
        if image is None:
            image = self._image_generator.create_empty_map_image()
        _LOGGER.debug("img: number of rooms: %d, numbers: %s",
                      len(rooms_raw), rooms_raw.keys())
        rooms = {}
        for number, room in rooms_raw.items():
            rooms[number] = Room(
                self.coord_transformer.image_to_map_x(room[0] + image_left),
                self.coord_transformer.image_to_map_y(room[1] + image_top),
                self.coord_transformer.image_to_map_x(room[2] + image_left),
                self.coord_transformer.image_to_map_y(room[3] + image_top),
                number,
            )
        # ImageData re-derives its own trim/scale from image_config + the
        # raw (untrimmed) width/height passed below, so the auto-detected
        # trim and scale computed by the image parser have to be threaded
        # through here too, otherwise path/wall/vacuum-position drawing
        # would be positioned/scaled against the configured values instead
        # of the ones actually used to size the rendered image.
        effective_image_config = replace(self._image_config, trim=effective_trim, scale=effective_scale)
        return (
            ImageData(
                image_size,
                image_top,
                image_left,
                image_height,
                image_width,
                effective_image_config,
                image,
                self.coord_transformer.map_to_image,
            ),
            rooms,
            cleaned_areas,
        )

    @staticmethod
    def _parse_history() -> Path:
        path_points = []
        for pt in IjaiMapDataParser.robot_map.historyPose.points:
            # 0: taxi, 1: working
            path_points.append(Point(x=pt.x, y=pt.y))
        return Path(len(path_points), 1, 0, [path_points])

    @staticmethod
    def _parse_restricted_areas() -> tuple[list[Wall], list[Area], list[Area]]:
        walls = []
        no_go_areas = []
        no_mop_areas = []

        for virtualWall in IjaiMapDataParser.robot_map.virtualWalls:
            p1, p2, p3, p4 = virtualWall.points

            if virtualWall.type == IjaiMapDataParser.VIRTUALWALL_TYPE_WALL:
                walls.append(Wall(p1.x, p1.y, p3.x, p3.y))
            elif virtualWall.type == IjaiMapDataParser.VIRTUALWALL_TYPE_NO_GO:
                no_go_areas.append(
                    Area(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, p4.x, p4.y))
            elif virtualWall.type == IjaiMapDataParser.VIRTUALWALL_TYPE_NO_MOP:
                no_mop_areas.append(
                    Area(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, p4.x, p4.y))
        return walls, no_go_areas, no_mop_areas

    @staticmethod
    def _parse_cleaning_zones() -> list[Zone]:
        zones = []
        for areaInfo in IjaiMapDataParser.robot_map.areasInfo:
            zones.append(Zone(areaInfo.points[0].x,
                              areaInfo.points[0].y,
                              areaInfo.points[2].x,
                              areaInfo.points[2].y))
        return zones

    @staticmethod
    def _parse_goto_point() -> Point | None:
        for navigationPoint in IjaiMapDataParser.robot_map.navigationPoints:
            if (
                    navigationPoint.status == 0
                    and navigationPoint.pointType == 1
                    and navigationPoint.x != 1100.0  # outside map
                    and navigationPoint.y != 1100.0):
                return Point(navigationPoint.x,
                             navigationPoint.y,
                             navigationPoint.phi * 180 / math.pi)
        return None

    @staticmethod
    def _parse_rooms(map_data_rooms: dict[int, Room]) -> None:
        map_id = IjaiMapDataParser.robot_map.mapHead.mapHeadId
        # While the vacuum is returning to the dock or idle, mapInfo may not
        # list the current map; the room data below is still valid.
        current_map = None
        for map_data in IjaiMapDataParser.robot_map.mapInfo:
            if map_data.mapHeadId == map_id:
                current_map = map_data
                break
        if current_map is not None:
            _LOGGER.debug("map#%d: %s", current_map.mapHeadId, current_map.mapName)
        for r in IjaiMapDataParser.robot_map.roomDataInfo:
            if map_data_rooms is not None and r.roomId in map_data_rooms:
                map_data_rooms[r.roomId].name = r.roomName
                map_data_rooms[r.roomId].pos_x = r.roomNamePost.x
                map_data_rooms[r.roomId].pos_y = r.roomNamePost.y

            room_text_pos = Point(r.roomNamePost.x, r.roomNamePost.y)
            _LOGGER.debug("room#%d: %s %s", r.roomId,
                          r.roomName, room_text_pos)
