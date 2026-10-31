"""Ijai map image parser."""

import logging

from PIL import Image
from PIL.Image import Image as ImageType
from PIL.Image import Resampling
from vacuum_map_parser_base.config.color import ColorsPalette, SupportedColor
from vacuum_map_parser_base.config.drawable import Drawable
from vacuum_map_parser_base.config.image_config import ImageConfig, TrimConfig
from vacuum_map_parser_base.map_data import Point

_LOGGER = logging.getLogger(__name__)

# Ijai vacuums report a fixed-size map canvas (e.g. 800x800) regardless of
# how much of it has actually been explored, so the explored area can end
# up as a small island surrounded by a large unexplored (MAP_OUTSIDE)
# background. Pad the auto-detected content bounding box by this many
# source pixels on each side before trimming to it.
AUTO_TRIM_PADDING = 15

# After auto-cropping, the explored area is often tiny in absolute pixel
# terms. If scaling it up to at least this many pixels on its longer edge
# would need more than the configured scale, bump the scale to reach it,
# so maps aren't blurry/pixelated by default.
AUTO_SCALE_TARGET_DIMENSION = 960


class IjaiImageParser:

    """Ijai map image parser."""

    MAP_OUTSIDE = 0x00
    MAP_WALL = 0xFF
    MAP_SCAN = 0x01
    MAP_NEW_DISCOVERED_AREA = 0x02
    MAP_ROOM_MIN = 10
    MAP_ROOM_MAX = 59
    MAP_SELECTED_ROOM_MIN = 60
    MAP_SELECTED_ROOM_MAX = 109

    def __init__(self, palette: ColorsPalette, image_config: ImageConfig, drawables: list[Drawable]):
        self._palette = palette
        self._image_config = image_config
        self._drawables = drawables
        self.color_map = {
            # Fully transparent rather than a filled color (matches Dreame's
            # rendering): the unexplored background should show through to
            # whatever is behind the image (card/page background), not
            # paint over it with a fixed color.
            IjaiImageParser.MAP_OUTSIDE: (255, 255, 255, 0),
            IjaiImageParser.MAP_WALL: palette.get_color(SupportedColor.MAP_WALL_V2),
            IjaiImageParser.MAP_SCAN: palette.get_color(SupportedColor.SCAN),
            IjaiImageParser.MAP_NEW_DISCOVERED_AREA: palette.get_color(SupportedColor.NEW_DISCOVERED_AREA)}

    def _compute_effective_trim(self, map_data: bytes, width: int, height: int) -> TrimConfig:
        configured = self._image_config.trim
        min_x = min_y = None
        max_x = max_y = -1
        for i, pixel_type in enumerate(map_data):
            if pixel_type != IjaiImageParser.MAP_OUTSIDE:
                x = i % width
                y = i // width
                if min_x is None or x < min_x:
                    min_x = x
                if min_y is None or y < min_y:
                    min_y = y
                max_x = max(max_x, x)
                max_y = max(max_y, y)

        if min_x is None or min_y is None:
            # Nothing explored yet; leave the configured trim as-is.
            return configured

        auto_left = max(0, min_x - AUTO_TRIM_PADDING) / width * 100
        auto_right = max(0, width - 1 - max_x - AUTO_TRIM_PADDING) / width * 100
        auto_bottom = max(0, min_y - AUTO_TRIM_PADDING) / height * 100
        auto_top = max(0, height - 1 - max_y - AUTO_TRIM_PADDING) / height * 100

        return TrimConfig(
            left=max(configured.left, auto_left),
            right=max(configured.right, auto_right),
            top=max(configured.top, auto_top),
            bottom=max(configured.bottom, auto_bottom),
        )

    def parse(
        self, map_data: bytes, width: int, height: int
    ) -> tuple[ImageType | None, dict[int, tuple[int, int, int, int]], set[int], TrimConfig, float]:
        rooms: dict[int, tuple[int, int, int, int]] = {}
        cleaned_areas = set()
        _LOGGER.debug("ijai parser: image_config = %s", self._image_config)
        configured_scale = self._image_config.scale
        effective_trim = self._compute_effective_trim(map_data, width, height)
        trim_left = int(effective_trim.left * width / 100)
        trim_right = int(effective_trim.right * width / 100)
        trim_top = int(effective_trim.top * height / 100)
        trim_bottom = int(effective_trim.bottom * height / 100)
        trimmed_height = height - trim_top - trim_bottom
        trimmed_width = width - trim_left - trim_right
        if trimmed_width == 0 or trimmed_height == 0:
            return None, {}, set(), effective_trim, configured_scale

        longer_edge = max(trimmed_width, trimmed_height)
        auto_scale = AUTO_SCALE_TARGET_DIMENSION / longer_edge if longer_edge > 0 else 1
        scale = max(configured_scale, auto_scale)
        _LOGGER.debug(
            "trimmed_width=%d trimmed_height=%d configured_scale=%s auto_scale=%.3f "
            "effective_scale=%.3f final_image_size=%dx%d",
            trimmed_width, trimmed_height, configured_scale, auto_scale, scale,
            int(trimmed_width * scale), int(trimmed_height * scale),
        )

        image = Image.new('RGBA', (trimmed_width, trimmed_height))
        pixels = image.load()
        _LOGGER.debug("trim_bottom = %s, trim_top = %s, trim_left = %s, trim_right = %s",
                      trim_bottom, trim_top, trim_left, trim_right)
        unknown_pixels = set()
        for img_y in range(trimmed_height):
            y = trimmed_height - 1 - img_y
            for img_x in range(trimmed_width):
                x = img_x
                pixel_type = map_data[(img_y + trim_bottom)
                                      * width + x + trim_left]
                if pixel_type in self.color_map:
                    pixels[x, y] = self.color_map[pixel_type]
                elif IjaiImageParser.MAP_ROOM_MIN <= pixel_type <= IjaiImageParser.MAP_SELECTED_ROOM_MAX:
                    room_x = img_x + trim_left
                    room_y = img_y + trim_bottom
                    room_number = pixel_type
                    if pixel_type >= IjaiImageParser.MAP_SELECTED_ROOM_MIN:
                        room_number = pixel_type - IjaiImageParser.MAP_SELECTED_ROOM_MIN + \
                            IjaiImageParser.MAP_ROOM_MIN
                        cleaned_areas.add(room_number)
                    rooms[room_number] = (room_x, room_y, room_x, room_y) \
                        if room_number not in rooms \
                        else (min(rooms[room_number][0], room_x),
                              min(rooms[room_number][1], room_y),
                              max(rooms[room_number][2], room_x),
                              max(rooms[room_number][3], room_y))
                    pixels[x, y] = self._palette.get_room_color(room_number)
                else:
                    pixels[x, y] = self._palette.get_color(
                        SupportedColor.UNKNOWN)
                    unknown_pixels.add(pixel_type)
                    _LOGGER.debug(
                        "unknown pixel [%s,%s] = %s", x, y, pixel_type)
        if scale != 1 and trimmed_width != 0 and trimmed_height != 0:
            # LANCZOS blurs the hard edges between rooms/walls on this
            # flat-color map; NEAREST keeps them crisp. The bigger win was
            # the auto-crop putting more of the image where the floor plan
            # actually is, not smoothing the resize.
            image = image.resize(
                (int(trimmed_width * scale), int(trimmed_height * scale)), resample=Resampling.NEAREST)
        if len(unknown_pixels) > 0:
            _LOGGER.warning('unknown pixel_types: %s', unknown_pixels)
        return image, rooms, cleaned_areas, effective_trim, scale

    @staticmethod
    def get_current_vacuum_room(map_data: bytes, vacuum_position_on_image: Point, image_width: int) -> int | None:
        _LOGGER.debug("pos on image: %s", vacuum_position_on_image)
        pixel_type = map_data[int(vacuum_position_on_image.y)
                              * image_width + int(vacuum_position_on_image.x)]
        if IjaiImageParser.MAP_ROOM_MIN <= pixel_type <= IjaiImageParser.MAP_ROOM_MAX:
            return pixel_type
        if IjaiImageParser.MAP_SELECTED_ROOM_MIN <= pixel_type <= IjaiImageParser.MAP_SELECTED_ROOM_MAX:
            return pixel_type - IjaiImageParser.MAP_SELECTED_ROOM_MIN + IjaiImageParser.MAP_ROOM_MIN
        return None
