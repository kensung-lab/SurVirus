from random import randint

class RandomPositionGenerator:
    def __init__(self, sampling_regions):
        self.sampling_regions = sampling_regions
        self.reference_len = sum([r[2]-r[1] for r in sampling_regions])

    def next(self):
        # according to documentation this is supposed to generate within [a, b), but sometimes it generates b for me
        r = randint(0, self.reference_len-1)
        for region in self.sampling_regions:
            k, start, stop = region
            if r < stop-start:
                chr = k
                break
            r -= stop-start
        return (chr, r+start)
