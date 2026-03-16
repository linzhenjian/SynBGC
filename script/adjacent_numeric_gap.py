# adjacent_numeric_gap.py
import re
import statistics

def adjacent_numeric_gap(genes):
    """Return the median gap size between adjacent genes."""
    nums = []
    for g in genes:
        match = re.search(r'_g(\d+)$', g)
        if match:
            nums.append(int(match.group(1)))

    gaps = []
    for i in range(len(nums) - 1):
        diff = nums[i+1] - nums[i]
        if diff > 1:
            gaps.append(diff)

    if gaps:
        return statistics.median(gaps)
    else:
        return None

