import pandas as pd
df = pd.read_csv('./J010534p234960.coarse.errfix.bin6.mwcorr.ascii',delim_whitespace=True)

lya_emitter_subset = df[df['EW'] >= 20]

