"""Ijai map parser."""

import logging
import math
import zlib
from typing import Any

from vacuum_map_parser_base.config.color import ColorsPalette
from vacuum_map_parser_base.config.drawable import Drawable
from vacuum_map_parser_base.config.image_config import ImageConfig
from vacuum_map_parser_base.config.size import Sizes
from vacuum_map_parser_base.config.text import Text
from vacuum_map_parser_base.map_data import Area, ImageData, MapData, Path, Point, Room, Wall, Zone
from vacuum_map_parser_base.map_data_parser import MapDataParser
import vacuum_map_parser_ijai.RobotMap_pb2 as RobotMap
import vacuum_map_parser_ijai.beautify_min as Beautify

from .image_parser import IjaiImageParser
from .aes_decryptor import decrypt

_LOGGER = logging.getLogger(__name__)


class IjaiMapDataParser(MapDataParser):
    """Ijai map parser."""

    POSITION_UNKNOWN = 1100
    VIRTUALWALL_TYPE_WALL = 2
    VIRTUALWALL_TYPE_NO_MOP = 6
    VIRTUALWALL_TYPE_NO_GO = 3

    robot_map = RobotMap.RobotMap()
    _map_to_img_scale = Point(0,0)

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

        IjaiMapDataParser.robot_map.ParseFromString(raw)
        map_head = IjaiMapDataParser.robot_map.mapHead
        IjaiMapDataParser._map_to_img_scale = Point(map_head.sizeX/(map_head.maxX - map_head.minX), 
                                                    map_head.sizeY/(map_head.maxY - map_head.minY))

        if hasattr(self.robot_map, "mapData"):
            map_data.image, map_data.rooms, map_data.cleaned_rooms = self._parse_image()

        if hasattr(self.robot_map, "historyPose"):
            map_data.path = IjaiMapDataParser._parse_history()

        if hasattr(self.robot_map, "chargeStation"):
            pos_info = self.robot_map.chargeStation
            map_data.charger = Point(x = pos_info.x, y = pos_info.y, a = pos_info.phi * 180 / math.pi)
            _LOGGER.debug("pos: %s", map_data.charger)

        if hasattr(self.robot_map, "currentPose"):
            pos_info = self.robot_map.currentPose
            map_data.vacuum_position = Point(x = pos_info.x, y = pos_info.y, a = pos_info.phi * 180 / math.pi)
            _LOGGER.debug("pos: %s", map_data.vacuum_position)

        if hasattr(self.robot_map, "mapInfo") and hasattr(self.robot_map, "roomDataInfo") and map_data.rooms is not None:
            IjaiMapDataParser._parse_rooms(map_data.rooms)

        if hasattr(self.robot_map, "virtualWalls"):
            map_data.walls, map_data.no_go_areas, map_data.no_mopping_areas = IjaiMapDataParser._parse_restricted_areas()

        if map_data.rooms is not None:
            _LOGGER.debug("rooms: %s", [str(room) for number, room in map_data.rooms.items()])
            if map_data.rooms is not None and len(map_data.rooms) > 0 and map_data.vacuum_position is not None:
                vacuum_position_on_image = IjaiMapDataParser._map_to_image(map_data.vacuum_position)
                map_data.vacuum_room = IjaiImageParser.get_current_vacuum_room(self.robot_map.mapData.mapData, vacuum_position_on_image, IjaiMapDataParser.robot_map.mapHead.sizeX)
                _LOGGER.debug(f"MinX={IjaiMapDataParser.robot_map.mapHead.minX} MaxX={IjaiMapDataParser.robot_map.mapHead.maxX}")
                if map_data.vacuum_room is not None:
                    map_data.vacuum_room_name = map_data.rooms[map_data.vacuum_room].name
                _LOGGER.debug("current vacuum room: %s", map_data.vacuum_room)
        return map_data

    @staticmethod
    def _map_to_image(p: Point) -> Point:
        scale = 1 / IjaiMapDataParser.robot_map.mapHead.resolution
        return Point((p.x - IjaiMapDataParser.robot_map.mapHead.minX) * IjaiMapDataParser._map_to_img_scale.x, 
                     (p.y - IjaiMapDataParser.robot_map.mapHead.minY) * IjaiMapDataParser._map_to_img_scale.y)

    @staticmethod
    def _image_to_map(x: float) -> float:
        return (x/IjaiMapDataParser._map_to_img_scale.x + IjaiMapDataParser.robot_map.mapHead.minX) 

    def _parse_image(self) -> tuple[ImageData, dict[int, Room], set[int]]:
        image_left = 0
        image_top = 0
        image_width = self.robot_map.mapHead.sizeX
        image_height = self.robot_map.mapHead.sizeY
        image_size = image_height * image_width
        _LOGGER.debug("width: %d, height: %d", image_width, image_height)

        mapData_temp = self.robot_map.mapData.mapData

        if (len(set(mapData_temp).symmetric_difference([0, 128, 127])) == 0 and len(self.robot_map.roomChain) > 0 and self.robot_map.mapType == 0):
            buautify_obj = Beautify.BeautifyMap()
            buautify_obj.setMap(self.robot_map)
            buautify_obj.transform()
            buautify_obj.roomColorByChain(self.robot_map.roomChain)
            buautify_obj.fillInternalObstacles()

            mapData_temp = buautify_obj.getMap()

            for i in range(len(mapData_temp)):
                if mapData_temp[i] < 0:
                    mapData_temp[i] = (256 + mapData_temp[i]) % 256
                elif mapData_temp[i] > 255:
                    mapData_temp[i] = mapData_temp[i] % 256
                elif mapData_temp[i] == 30:
                    mapData_temp[i] = 0
                elif mapData_temp[i] == 40:
                    mapData_temp[i] = 255
            self.robot_map.mapData.mapData = bytes(mapData_temp)

        image, rooms_raw, cleaned_areas, cleaned_areas_layer = self._image_parser.parse(self.robot_map.mapData.mapData, image_width, image_height)
        if image is None:
            image = self._image_generator.create_empty_map_image()
        _LOGGER.debug("img: number of rooms: %d, numbers: %s", len(rooms_raw), rooms_raw.keys())
        rooms = {}
        for number, room in rooms_raw.items():
            rooms[number] = Room(
                IjaiMapDataParser._image_to_map(room[0] + image_left),
                IjaiMapDataParser._image_to_map(room[1] + image_top),
                IjaiMapDataParser._image_to_map(room[2] + image_left),
                IjaiMapDataParser._image_to_map(room[3] + image_top),
                number,
            )
        return (
            ImageData(
                image_size,
                image_top,
                image_left,
                image_height,
                image_width,
                self._image_config,
                image,
                IjaiMapDataParser._map_to_image,
                additional_layers={Drawable.CLEANED_AREA: cleaned_areas_layer},
            ),
            rooms,
            cleaned_areas,
        )

    @staticmethod
    def _parse_history() -> Path:
        path_points = []
        for pt in IjaiMapDataParser.robot_map.historyPose.points:
            # 0: taxi, 1: working
            path_points.append(Point(x = pt.x, y = pt.y))
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
    def _parse_rooms(map_data_rooms: dict[int, Room]) -> None:
        map_id = IjaiMapDataParser.robot_map.mapHead.mapHeadId
        for map_data in IjaiMapDataParser.robot_map.mapInfo:
            if (map_data.mapHeadId == map_id):
                current_map = map_data
                break
        map_name = current_map.mapName
        _LOGGER.debug("map#%d: %s", current_map.mapHeadId, map_name)
        for r in IjaiMapDataParser.robot_map.roomDataInfo:
            if map_data_rooms is not None and r.roomId in map_data_rooms:
                map_data_rooms[r.roomId].name = r.roomName
                map_data_rooms[r.roomId].pos_x = r.roomNamePost.x
                map_data_rooms[r.roomId].pos_y = r.roomNamePost.y

            room_text_pos = Point(r.roomNamePost.x, r.roomNamePost.y)
            _LOGGER.debug("room#%d: %s %s", r.roomId, r.roomName, room_text_pos)