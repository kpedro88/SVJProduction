flist = []
for width in [0, 0.01, 0.05, 0.10, 0.15, 0.20, 0.30]:
    flist.append({
        "channel": "s",
        "mMediator": 3100,
        "mDark": 20,
        "rinv": 0.3,
        "alpha": "peak",
        "width": width,
    })
