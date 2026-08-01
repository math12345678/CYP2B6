"""
Sanity checks that the H-bond data is real MD output, not fabricated/mock
values. Real trajectory data should show frame-to-frame fluctuation (not a
constant value), a plausible distribution shape, and consistency with the
independently-known frame count (30001) and mutation-site chemistry already
confirmed in verify_hbond_numbering.py.

Run from ~/Desktop/Research/Research_Projects/CYP2B6, inside the cyp2b6 env.
"""
import numpy as np

def load_count(path):
    return np.loadtxt(path, delimiter=",", skiprows=1)

print("=== Per-frame H-bond count statistics (WT) ===")
data = load_count("WT/hbond_count_WT.csv")
frames, counts = data[:, 0], data[:, 1]
print(f"n_frames: {len(frames)} (expect 30001, matching the trajectory length reported earlier)")
print(f"frame range: {frames.min():.0f} to {frames.max():.0f}")
print(f"count mean: {counts.mean():.2f}, std: {counts.std():.2f}, min: {counts.min():.0f}, max: {counts.max():.0f}")
print(f"unique count values seen: {len(set(counts))} (a constant/mock series would show 1)")

print("\n=== First 10 and last 10 frame counts (WT) -- should show real fluctuation, not a flat line ===")
print("first 10:", counts[:10].astype(int).tolist())
print("last 10:", counts[-10:].astype(int).tolist())

print("\n=== Autocorrelation sanity check ===")
# Real MD data is autocorrelated frame-to-frame (physical continuity), unlike
# purely random/fabricated noise. Check correlation between count[t] and
# count[t+1].
c0 = counts[:-1]
c1 = counts[1:]
corr = np.corrcoef(c0, c1)[0, 1]
print(f"lag-1 autocorrelation: {corr:.3f} (real MD trajectories are typically > 0.3-0.5; "
      f"independent random noise would be close to 0)")

print("\n=== Cross-check against a second, independent system (P428T) ===")
data2 = load_count("P428T/hbond_count_P428T.csv")
counts2 = data2[:, 1]
print(f"P428T count mean: {counts2.mean():.2f}, std: {counts2.std():.2f} "
      f"(should be in the same ballpark as WT's {counts.mean():.2f}, but not identical)")
print(f"WT and P428T count arrays identical?: {np.array_equal(counts, counts2)} "
      f"(should be False -- identical arrays across different systems would indicate copy-paste/mock data)")
