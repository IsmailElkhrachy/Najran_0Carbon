"""Minimal ESRI shapefile writer (Polygon, type 5) with DBF and PRJ."""
import struct
PRJ_UTM38N=('PROJCS["WGS_1984_UTM_Zone_38N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",'
 'SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
 'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],'
 'PARAMETER["Central_Meridian",45.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')
def _bbox(parts):
    xs=[x for p in parts for x,_ in p]; ys=[y for p in parts for _,y in p]; return min(xs),min(ys),max(xs),max(ys)
def write_polygons(base, polys, fields, records, prj=PRJ_UTM38N):
    """polys: list of parts (list of rings [(x,y),...], outer CW, holes CCW, closed).
    fields: list of (name, 'N'|'C'|'F', length, decimals). records: list of dicts."""
    recs=[]
    for parts in polys:
        npts=sum(len(p) for p in parts); bb=_bbox(parts)
        idx=[];c=0
        for p in parts: idx.append(c); c+=len(p)
        body=struct.pack('<i4d2i',5,*bb,len(parts),npts)+struct.pack(f'<{len(parts)}i',*idx)
        body+=b''.join(struct.pack('<2d',x,y) for p in parts for x,y in p)
        recs.append(body)
    allbb=[_bbox(p) for p in polys]; BB=(min(b[0] for b in allbb),min(b[1] for b in allbb),max(b[2] for b in allbb),max(b[3] for b in allbb))
    def header(flen): return struct.pack('>7i',9994,0,0,0,0,0,flen)+struct.pack('<2i4d4d',1000,5,*BB,0,0,0,0)
    shp=b''.join(struct.pack('>2i',i+1,len(r)//2)+r for i,r in enumerate(recs))
    with open(base+'.shp','wb') as f: f.write(header((100+len(shp))//2)); f.write(shp)
    off=50; shx=b''
    for r in recs: shx+=struct.pack('>2i',off,len(r)//2); off+=4+len(r)//2
    with open(base+'.shx','wb') as f: f.write(header((100+len(shx))//2)); f.write(shx)
    rl=1+sum(f[2] for f in fields)
    hdr=struct.pack('<B3BIHH20x',3,126,1,1,len(records),32+32*len(fields)+1,rl)
    fd=b''.join(struct.pack('<11sc4xBB14x',n.encode()[:11],t.encode(),l,d) for n,t,l,d in fields)
    rows=b''
    for r in records:
        row=b' '
        for n,t,l,d in fields:
            v=r.get(n,'')
            if t in 'NF': s=(f'{v:.{d}f}' if v is not None else '').rjust(l)
            else: s=str(v)
            b=s.encode('utf-8')[:l]; row+=b.ljust(l) if t=='C' else b.rjust(l)
        rows+=row
    with open(base+'.dbf','wb') as f: f.write(hdr+fd+b'\r'+rows+b'\x1a')
    open(base+'.prj','w').write(prj); open(base+'.cpg','w').write('UTF-8')
