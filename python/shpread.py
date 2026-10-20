"""Minimal pure-Python reader for ESRI shapefiles (.shp + .dbf), plus
UTM->WGS84 conversion. Used because GeoPandas is not available here."""
import math
import struct


def read_dbf(path, encoding="utf-8"):
    with open(path, "rb") as f:
        hdr = f.read(32)
        nrec, hlen, rlen = struct.unpack("<IHH", hdr[4:12])
        fields = []
        while True:
            fd = f.read(32)
            if fd[0] == 0x0D:
                break
            name = fd[:11].split(b"\x00")[0].decode("ascii", "replace")
            ftype = chr(fd[11])
            flen = fd[16]
            fdec = fd[17]
            fields.append((name, ftype, flen, fdec))
        f.seek(hlen)
        rows = []
        for _ in range(nrec):
            rec = f.read(rlen)
            if not rec:
                break
            pos = 1
            row = {}
            for name, ftype, flen, fdec in fields:
                raw = rec[pos:pos + flen]
                pos += flen
                txt = raw.decode(encoding, "replace").strip()
                if ftype in "NF":
                    try:
                        row[name] = float(txt) if txt else None
                    except ValueError:
                        row[name] = None
                else:
                    row[name] = txt
            rows.append(row)
    return fields, rows


def read_shp(path):
    """Return (shape_type, list of geometries). Each geometry is a list of
    parts; each part a list of (x, y)."""
    geoms = []
    with open(path, "rb") as f:
        hdr = f.read(100)
        stype = struct.unpack("<i", hdr[32:36])[0]
        while True:
            rh = f.read(8)
            if len(rh) < 8:
                break
            _, clen = struct.unpack(">ii", rh)
            content = f.read(clen * 2)
            t = struct.unpack("<i", content[:4])[0]
            if t == 0:
                geoms.append([])
                continue
            if t in (1, 11, 21):
                x, y = struct.unpack("<2d", content[4:20])
                geoms.append([[(x, y)]])
                continue
            nparts, npts = struct.unpack("<2i", content[36:44])
            parts = list(struct.unpack(f"<{nparts}i", content[44:44 + 4 * nparts]))
            off = 44 + 4 * nparts
            pts = struct.unpack(f"<{2 * npts}d", content[off:off + 16 * npts])
            xy = list(zip(pts[0::2], pts[1::2]))
            parts.append(npts)
            geoms.append([xy[parts[i]:parts[i + 1]] for i in range(nparts)])
    return stype, geoms


def ring_area(ring):
    a = 0.0
    for (x1, y1), (x2, y2) in zip(ring, ring[1:] + ring[:1]):
        a += x1 * y2 - x2 * y1
    return a / 2.0


def polygon_area(parts):
    """ESRI rings: outer clockwise (negative shoelace), holes counter-clockwise."""
    return -sum(ring_area(r) for r in parts)


def line_length(parts):
    return sum(math.hypot(x2 - x1, y2 - y1)
               for p in parts for (x1, y1), (x2, y2) in zip(p, p[1:]))


def utm_to_wgs84(e, n, zone=38, north=True):
    a = 6378137.0
    f = 1 / 298.257223563
    k0 = 0.9996
    e2 = f * (2 - f)
    ep2 = e2 / (1 - e2)
    x = e - 500000.0
    y = n if north else n - 10000000.0
    m = y / k0
    mu = m / (a * (1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256))
    e1 = (1 - math.sqrt(1 - e2)) / (1 + math.sqrt(1 - e2))
    p1 = (mu + (3 * e1 / 2 - 27 * e1 ** 3 / 32) * math.sin(2 * mu)
          + (21 * e1 ** 2 / 16 - 55 * e1 ** 4 / 32) * math.sin(4 * mu)
          + (151 * e1 ** 3 / 96) * math.sin(6 * mu))
    n1 = a / math.sqrt(1 - e2 * math.sin(p1) ** 2)
    t1 = math.tan(p1) ** 2
    c1 = ep2 * math.cos(p1) ** 2
    r1 = a * (1 - e2) / (1 - e2 * math.sin(p1) ** 2) ** 1.5
    d = x / (n1 * k0)
    lat = p1 - (n1 * math.tan(p1) / r1) * (
        d ** 2 / 2 - (5 + 3 * t1 + 10 * c1 - 4 * c1 ** 2 - 9 * ep2) * d ** 4 / 24
        + (61 + 90 * t1 + 298 * c1 + 45 * t1 ** 2 - 252 * ep2 - 3 * c1 ** 2) * d ** 6 / 720)
    lon = (d - (1 + 2 * t1 + c1) * d ** 3 / 6
           + (5 - 2 * c1 + 28 * t1 - 3 * c1 ** 2 + 8 * ep2 + 24 * t1 ** 2) * d ** 5 / 120) / math.cos(p1)
    return math.degrees(lon) + (zone * 6 - 183), math.degrees(lat)


def wgs84_to_utm(lon, lat, zone=38):
    a = 6378137.0
    f = 1 / 298.257223563
    k0 = 0.9996
    e2 = f * (2 - f)
    ep2 = e2 / (1 - e2)
    lat_r = math.radians(lat)
    lon0 = math.radians(zone * 6 - 183)
    n = a / math.sqrt(1 - e2 * math.sin(lat_r) ** 2)
    t = math.tan(lat_r) ** 2
    c = ep2 * math.cos(lat_r) ** 2
    aa = math.cos(lat_r) * (math.radians(lon) - lon0)
    m = a * ((1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256) * lat_r
             - (3 * e2 / 8 + 3 * e2 ** 2 / 32 + 45 * e2 ** 3 / 1024) * math.sin(2 * lat_r)
             + (15 * e2 ** 2 / 256 + 45 * e2 ** 3 / 1024) * math.sin(4 * lat_r)
             - (35 * e2 ** 3 / 3072) * math.sin(6 * lat_r))
    e = k0 * n * (aa + (1 - t + c) * aa ** 3 / 6
                  + (5 - 18 * t + t ** 2 + 72 * c - 58 * ep2) * aa ** 5 / 120) + 500000.0
    nn = k0 * (m + n * math.tan(lat_r) * (aa ** 2 / 2 + (5 - t + 9 * c + 4 * c ** 2) * aa ** 4 / 24
                                          + (61 - 58 * t + t ** 2 + 600 * c - 330 * ep2) * aa ** 6 / 720))
    return e, nn


def point_in_ring(x, y, ring):
    inside = False
    j = len(ring) - 1
    for i in range(len(ring)):
        xi, yi = ring[i]
        xj, yj = ring[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi + 1e-300) + xi):
            inside = not inside
        j = i
    return inside


def point_in_polygon(x, y, parts):
    return sum(point_in_ring(x, y, r) for r in parts) % 2 == 1
