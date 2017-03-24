function y = uint32todouble(x, minx, maxx, miny, maxy)
y = x*(maxy-miny)/double(maxx - minx);