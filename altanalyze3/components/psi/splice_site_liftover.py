"""Auditable point liftover of ordered, one-based splice-site pairs.

Uses pyliftover's native point conversion. All chain alternatives remain visible;
the two endpoints of an accepted junction must share an alignment chain.
"""
from functools import lru_cache
import re


def parse_coordinate_junction(text):
    match = re.fullmatch(r'(chr[^:]+):(\d+)-(\d+)', str(text))
    if not match:
        raise ValueError(f'Unresolved coordinate junction: {text}')
    chrom, first, second = match.groups()
    first, second = int(first), int(second)
    if min(first, second) < 1:
        raise ValueError('Expected one-based genomic nucleotide positions')
    return chrom, first, second


class SpliceSiteLiftover:
    def __init__(self, chain_file):
        from pyliftover import LiftOver
        self.native = LiftOver(str(chain_file))

    @lru_cache(maxsize=None)
    def point(self, chrom, coordinate_1based):
        """All native mappings, with their original alignment-chain identities."""
        position = int(coordinate_1based) - 1
        if position < 0:
            raise ValueError('Expected positive one-based coordinate')
        values = self.native.convert_coordinate(chrom, position, '+') or []
        blocks = self.native.chain_file.query(chrom, position) or []
        # convert_coordinate uses this identical stable score ordering; attach IDs
        # from the very same queried native blocks without reimplementing conversion.
        blocks = sorted(blocks, key=lambda b: b[2][1].score, reverse=True)
        assert len(blocks) == len(values)
        out = []
        for value, block in zip(values, blocks):
            target_chrom, target_zero, orientation, score = value
            chain = block[2][1]
            assert (target_chrom, orientation, score) == (chain.target_name, chain.target_strand, chain.score)
            out.append((target_chrom, int(target_zero) + 1, orientation, int(score), str(chain.id)))
        return tuple(out)

    def junction(self, text):
        chrom, first, second = parse_coordinate_junction(text)
        a, b = self.point(chrom, first), self.point(chrom, second)
        mapped = []
        for x in a:
            for y in b:
                if x[0] == y[0] and x[2:] == y[2:]:
                    mapped.append(dict(chromosome=x[0], first=x[1], second=y[1],
                        chain_orientation=x[2], chain_score=x[3], chain_id=x[4]))
        return dict(source_chromosome=chrom, source_first=first, source_second=second,
                    first_point_mappings=a, second_point_mappings=b,
                    same_chain_junction_mappings=mapped)
