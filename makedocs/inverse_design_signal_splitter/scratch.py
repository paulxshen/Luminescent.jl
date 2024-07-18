from cairosvg import svg2png
import textwrap
import shapely


def s2svg(area, bbox, dx):
    # specify margin in coordinate units

    width = bbox[2] - bbox[0]
    height = bbox[3] - bbox[1]
    props = {
        'version': '1.1',
        'baseProfile': 'full',
        'width': '{width:.0f}px'.format(width=round(width/dx)),
        'height': '{height:.0f}px'.format(height=round(height/dx)),
        'viewBox': '%.1f,%.1f,%.1f,%.1f' % (bbox[0], bbox[1], width, height),
        'xmlns': 'http://www.w3.org/2000/svg',
        'xmlns:ev': 'http://www.w3.org/2001/xml-events',
        'xmlns:xlink': 'http://www.w3.org/1999/xlink',
        # "fill": "#FFFFFF",
    }

    return textwrap.dedent(r'''
        <?xml version="1.0" encoding="utf-8" ?>
        <svg {attrs:s}>
        {data:s}
        </svg>
    ''').format(
        attrs=' '.join(['{key:s}="{val:s}"'.format(
            key=key, val=props[key]) for key in props]),
        data=area.svg()
    ).strip().replace('stroke-width="2.0"', 'stroke-width="0.0"')


p = shapely.Polygon([(0.1, 0.1), (1.1, 0.1), (1.1, 1.1), (0.1, 1.1)])
svg = s2svg(p, [0, 0, 2, 2], .2)
svg2png(bytestring=svg, write_to='_.png', background_color="#00000000",
        output_width=10, output_height=10)
