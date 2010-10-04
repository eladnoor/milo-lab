from groups import GroupContribution
from common import R
from pylab import log, zeros, vstack, pinv, dot, inf, plot, matrix, show, figure, NaN, isnan, find
import csv
import kegg

G = GroupContribution(sqlite_name="gibbs.sqlite", html_name="acetogens")
G.init()

reactions = [] # (RID, EC, sparse-reaction, dG0_r, pH, I, T

# Drake 2006
#reactions.append([134, '1.2.1.43', {11:-1, 5:-1, 58:1, 6:1}, 22, 7.0, 0, 300])
reactions.append([934, '6.3.4.3', {101:-1, 58:-1, 2:-1, 8:1, 9:1, 234:1}, -8, 7.0, 0, 300])
reactions.append([1655, '3.5.4.9', {234:-1, 445:1, 1:1}, -4, 7.0, 0, 300])
reactions.append([1220, '1.5.1.5', {445:-1, 5:-1, 143:1, 6:1}, -5, 7.0, 0, 300])
reactions.append([1224, '1.5.1.20', {143:-1, 5:-1, 440:1, 6:1}, -22, 7.0, 0, 300])

# NIST database
reactions.append([134, '1.2.1.43', {288:1, 5:1, 58:-1, 6:-1, 1:-1}, -R * 328 * log(650), 7.5, 0, 328]) # Yamamoto 1983, Buffer: triethanolamine-maleate (0.1 M)
reactions.append([934, '6.3.4.3', {101:-1, 58:-1, 2:-1, 8:1, 9:1, 234:1}, -R * 310 * log(41), 7.7, 0, 310]) # Himes 1962, Buffer: triethanolamine (0.1 M)
reactions.append([1655, '3.5.4.9', {234:1, 445:-1, 1:-1}, -R * 298 * log(11), 7.0, 0, 298]) # Kay 1960, Buffer: acetate
reactions.append([1655, '3.5.4.9', {234:1, 445:-1, 1:-1}, -R * 298 * log(1.84), 6.5, 0, 298]) # Lombrozo 1967, Buffer: potassium citrate (0.11 M)
reactions.append([1655, '3.5.4.9', {234:1, 445:-1, 1:-1}, -R * 298 * log(4.2), 6.5, 0, 298]) # Greenberg 1963, Buffer: potassium maleate (1.0 M) 
reactions.append([1655, '3.5.4.9', {234:1, 445:-1, 1:-1}, -R * 298 * log(50), 7.0, 0, 298]) # Suzuki 1973, Buffer: potassium maleate
reactions.append([1220, '1.5.1.5', {445:1, 5:1, 143:-1, 6:-1}, -R * 298 * log(0.14), 6.9, 0, 298]) # Uyeda 1967, Buffer: potassium maleate (0.05 M)
reactions.append([1220, '1.5.1.5', {445:1, 5:1, 143:-1, 6:-1}, -R * 303 * log(16), 7.3, 0, 303]) # Pelletier 1995, Buffer: potassium phosphate (0.025 M)
reactions.append([8550, '1.8.1.4', {2972:-1, 3:-1, 2051:1, 4:1}, -R * 295 * log(0.21), 7.1, 0, 295]) # Sanadi 1957, Buffer: phosphate (0.026 M) 
reactions.append([8550, '1.8.1.4', {2972:-1, 3:-1, 2051:1, 4:1}, -R * 295 * log(0.13), 7.1, 0, 295]) # Sanadi 1959, Buffer: phosphate (0.05 M)
reactions.append([945, '2.1.2.1', {143:-1, 37:-1, 1:-1, 101:1, 65:1}, -R * 310 * log(0.067), 7.4, 0, 310]) # Besson 1993, Buffer: KH2PO4 (0.02 M)

# Directly taken from Liegel 1985 Thesis
reactions.append([9093, None, {101:-1, 67:-1, 143:1}, -R * 311 * log(3e4), 7.0, 0.25, 311])
reactions.append([8550, '1.8.1.4', {2972:-1, 3:-1, 2051:1, 4:1}, -R * 298 * log(0.267), 7.08, 0.25, 298]) # Buffer: potassium phosphate (0.050 M) or sodium pyrophosphate (0.030 M) 
reactions.append([8550, '1.8.1.4', {2972:-1, 3:-1, 2051:1, 4:1}, -R * 311 * log(0.266), 7.06, 0.25, 311]) # Buffer: potassium phosphate (0.050 M) or sodium pyrophosphate (0.030 M)
reactions.append([3425, '1.4.4.2', {1242:1, 11:1, 37:-1, 2051:-1}, -R * 311 * log(0.03), 7.0, 0.25, 311])
#reactions.append([4125, '2.1.2.10', {1242:-1, 101:-1, 2972:1, 143:1, 14:1}, None, None, None, None])
reactions.append([1221, None, {143:1, 4:1, 14:1, 11:1, 101:-1, 3:-1, 37:-1}, -R * 311 * log(1.56e-3), 7.0, 0.25, 311]) # combined 2.1.2.10 + 1.4.4.2 + 1.8.1.4

# reactions for which we wish to predict the dG
reactions.append([4125, '2.1.2.10', {1242:1, 101:1, 2972:-1, 143:-1, 14:-1}, NaN, 7.0, 0.25, 311]) # DLP + 5,10-methylene-THF + NH3 => AMDLP + THF 
reactions.append([3425, '1.4.4.2', {1242:-1, 11:-1, 37:1, 2051:1}, NaN, None, None, None]) # AMDLP + CO2 => LP + glycine
reactions.append([8550, '1.8.1.4', {2972:1, 3:1, 2051:-1, 4:-1}, NaN, None, None, None]) # LP + NADH => DLP + NAD+
reactions.append([1221, None, {143:1, 4:1, 14:1, 11:1, 101:-1, 3:-1, 37:-1}, NaN, None, None, None]) # 5,10-methylene-THF + NH3 + CO2 + NADH => THF + glycine + NAD+
reactions.append([None, None, {440:-1, 10:-1, 11:-1, 5:-1, 101:1, 6:1, 24:1, 1:1}, NaN, 7.0, 0.25, 311]) # 5-methyl-THF + CoA + NADPH + CO2 -> THF + NADP+ + acetyl-CoA + H2O 
reactions.append([2289, '2.1.1.-', {440:-1, 6021:-1, 101:1, 6020:1}, NaN, None, None, None])
reactions.append([8433, '2.3.1.169', {237:-1, 6020:-1, 10:-1, 24:1, 6021:1}, NaN, None, None, None])
reactions.append([None, None, {58:-1, 2:-1, 4:-2, 11:-1, 14:-1, 37:1, 8:1, 9:1, 3:2, 1:1}, NaN, None, None, None])

all_cids = set()
for r in reactions:
    all_cids = all_cids.union(r[2].keys())
all_cids = sorted(list(all_cids))

S = zeros((len(reactions), len(all_cids)))
dG0_r = zeros((len(reactions), 1))
for r in xrange(len(reactions)):
    for (cid, coeff) in reactions[r][2].iteritems():
        c = all_cids.index(cid)
        S[r, c] = coeff
        dG0_r[r, 0] = reactions[r][3]

# move known formation energies to the other side
known_columns = []
dG0_f = zeros((len(all_cids), 1))
for c in xrange(len(all_cids)):
    if (all_cids[c] in [101, 2051]): # define as 0 (THF and LP)
        known_columns.append(c)
        dG0_f[c, 0] = 0
    elif (all_cids[c] in [440, 234, 143, 445, 2972, 1242, 6020, 6021]): # although GC can calculate it, I'd rather not use that value
        continue
    else:
        try:
            temp = G.estimate_dG0_keggcid(all_cids[c], pH=7.0, I=0.25, T=298.15)
            #temp = G.estimate_dG0_keggcid(all_cids[c], pH=7.0, I=0.0, T=311)
            known_columns.append(c)
            dG0_f[c, 0] = temp
        except Exception:
            continue

unknown_columns = sorted(list(set(range(len(all_cids))).difference(known_columns))) # find all the indices of columns not in known_columns
unknown_rows = find(isnan(dG0_r))
known_rows = sorted(list(set(range(len(reactions))).difference(unknown_rows))) # find all the indices of rows with measured dG0_r

S_measured = S[known_rows, :]
b = dG0_r[known_rows] - dot(S_measured[:, known_columns], dG0_f[known_columns])
S_red = S_measured[:, unknown_columns]

print "Formation energies from from GC and from linear regression: "
# linear regression to solve the missing formation energies (S*x = b)
inv_corr_mat = pinv(dot(S_red.T, S_red))
x = dot(dot(inv_corr_mat, S_red.T), b)
for c in xrange(len(unknown_columns)):
    dG0_f[unknown_columns[c], 0] = x[c, 0]

csv_out = csv.writer(open('../res/acetogens_compounds.csv', 'w'))
for c in xrange(len(all_cids)):
    if (c in known_columns):
        csv_out.writerow(['C%05d' % all_cids[c], '%8.1f' % dG0_f[c, 0], 'Group Contribution']) 
    else:
        csv_out.writerow(['C%05d' % all_cids[c], '%8.1f' % dG0_f[c, 0], 'Calculated'])

#figure()
#plot(b, dot(S_red, dG0_f[unknown_columns]), '.')
figure()
plot(dG0_r[known_rows], dot(S[known_rows, :], dG0_f), '.')
show()

csv_out = csv.writer(open('../res/acetogens_reactions.csv', 'w'))
csv_out.writerow(["RID", "EC", "REACTION", "dG0 measured", "dG0 calculated", "pH", "I", "T"])
for r in xrange(len(reactions)):
    (rid, ec, sparse, dG0, pH, I, T) = reactions[r]
    dG0_calculated = dot(S[r, :], dG0_f)
    if (r in known_rows):
        csv_out.writerow([rid, ec, kegg.unparse_reaction_formula(sparse), '%.1f' % dG0, '%.1f' % dG0_calculated, pH, I, T])
    else:
        csv_out.writerow([rid, ec, kegg.unparse_reaction_formula(sparse), 'N/A', '%.1f' % dG0_calculated, pH, I, T])

