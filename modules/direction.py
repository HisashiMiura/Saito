from enum import Enum


class Direction(Enum):

    S = 's'
    SW = 'sw'
    W = 'w'
    NW = 'nw'
    N = 'n'
    NE = 'ne'
    E = 'e'
    SE = 'se'
    TOP = 'top'
    BOTTOM = 'bottom'

    @property
    def alpha(self):

        match self:
            case Direction.S:
                return 0.0
            case Direction.SW:
                return 45.0
            case Direction.W:
                return 90.0
            case Direction.NW:
                return 135.0
            case Direction.N:
                return 180.0
            case Direction.NE:
                return 225.0
            case Direction.E:
                return 270.0
            case Direction.SE:
                return 315.0
            case Direction.TOP:
                return None
            case Direction.BOTTOM:
                return None
            case _:
                raise Exception()
    

