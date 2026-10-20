"""Read tiled, LZW-compressed GeoTIFFs (as exported by Earth Engine) without
imagecodecs: pure-Python TIFF-LZW decoder + tile assembly via tifffile."""
import numpy as np
import tifffile


def lzw_decode(data):
    out = bytearray()
    table = [bytes([i]) for i in range(256)] + [b"", b""]
    bitpos, nbits, prev = 0, 9, None
    total = len(data) * 8
    while bitpos + nbits <= total:
        byte = bitpos >> 3
        chunk = int.from_bytes(data[byte:byte + 4].ljust(4, b"\0"), "big")
        code = (chunk >> (32 - (bitpos & 7) - nbits)) & ((1 << nbits) - 1)
        bitpos += nbits
        if code == 256:                       # clear
            table = table[:258]
            nbits, prev = 9, None
            continue
        if code == 257:                       # end of information
            break
        if prev is None:
            entry = table[code]
        elif code < len(table):
            entry = table[code]
            table.append(prev + entry[:1])
        else:
            entry = prev + prev[:1]
            table.append(entry)
        out += entry
        prev = entry
        n = len(table)
        if n + 1 >= 2048:
            nbits = 12
        elif n + 1 >= 1024:
            nbits = 11
        elif n + 1 >= 512:
            nbits = 10
    return bytes(out)


def read_geotiff(path):
    """Return (array[bands, rows, cols], geotransform (x0, dx, y0, dy))."""
    with tifffile.TiffFile(path) as tf:
        p = tf.pages[0]
        tags = {t.name: t.value for t in p.tags.values()}
        h, w = p.imagelength, p.imagewidth
        spp = p.samplesperpixel
        dt = np.dtype(p.dtype)
        tw, th = p.tilewidth, p.tilelength
        fh = tf.filehandle
        out = np.zeros((spp, h, w), dtype=dt)
        tiles_x = (w + tw - 1) // tw
        tiles_y = (h + th - 1) // th
        planar = p.planarconfig == 2 or getattr(p.planarconfig, "value", 1) == 2
        for i, (off, cnt) in enumerate(zip(p.dataoffsets, p.databytecounts)):
            fh.seek(off)
            raw = lzw_decode(fh.read(cnt))
            pred = getattr(p.predictor, "value", p.predictor)
            if planar:
                band, ti = divmod(i, tiles_x * tiles_y)
                tile = np.frombuffer(raw[:tw * th * dt.itemsize], dtype=dt).reshape(th, tw).copy()
                if pred == 2:
                    tile = np.cumsum(tile, axis=1, dtype=dt)
                ty, tx = divmod(ti, tiles_x)
                y0, x0 = ty * th, tx * tw
                out[band, y0:y0 + th, x0:x0 + tw] = tile[:min(th, h - y0), :min(tw, w - x0)]
            else:
                tile = np.frombuffer(raw[:tw * th * spp * dt.itemsize], dtype=dt).reshape(th, tw, spp).copy()
                if pred == 2:
                    tile = np.cumsum(tile, axis=1, dtype=dt)
                ty, tx = divmod(i, tiles_x)
                y0, x0 = ty * th, tx * tw
                out[:, y0:y0 + th, x0:x0 + tw] = np.moveaxis(
                    tile[:min(th, h - y0), :min(tw, w - x0)], 2, 0)
        sc = tags["ModelPixelScaleTag"]
        tp = tags["ModelTiepointTag"]
        gt = (tp[3] - tp[0] * sc[0], sc[0], tp[4] + tp[1] * sc[1], -sc[1])
    return out, gt
