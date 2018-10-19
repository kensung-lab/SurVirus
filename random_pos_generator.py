from random import randint


class RandomPositionGenerator:
    def __init__(self, reference_fa):
        self.reference_fa = reference_fa
        self.reference_len = -1
        for k in reference_fa.keys():
            self.reference_len += len(reference_fa[k])

    def next(self):
        r = randint(0, self.reference_len)
        for k in self.reference_fa.keys():
            if r < len(self.reference_fa[k]):
                chr = k
                break
            r -= len(self.reference_fa[k])
        return (chr, r)
