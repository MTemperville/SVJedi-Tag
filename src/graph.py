class Orientation:
    FORWARD, REVERSE = "+", "-"

_FLIP = {"+": "-", "-": "+"}

class Graph:
    'Manual parsing of GFA file that computes orientation. Does not store sequences.'
    def __init__(self, gfa_file):
        self.segments = {}
        new = lambda: {"out": {"+": [], "-": []}}
        with open(gfa_file, encoding="UTF-8") as f:
            for line in f:
                t = line.rstrip("\n").split("\t", 5)
                if t[0] == "S":
                    self.segments.setdefault(t[1], new())
                elif t[0] == "L":
                    a, ao, b, bo = t[1:5]
                    self.segments.setdefault(a, new())["out"][ao].append((b, bo))
                    # reverse-complement edge, same as with_reverse_edges=True
                    self.segments.setdefault(b, new())["out"][_FLIP[bo]].append((a, _FLIP[ao]))