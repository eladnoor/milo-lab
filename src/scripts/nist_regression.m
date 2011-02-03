clear all;
S = load('../../res/nist/regress_S.txt');
dG0_r = load('../../res/nist/regress_dG0.txt');
cids = load('../../res/nist/regress_CID.txt');

nonzero_cols = find(sum(abs(S), 1)>0);
S = S(:, nonzero_cols);
cids = cids(nonzero_cols);
clear nonzero_cols;

dG0_f = S \ dG0_r;

figure(1);
plot(S * dG0_f, dG0_r, '.');

r_X = rref([S, dG0_r]);
r_S = r_X(:, 1:end-1);
r_dG0_r = r_X(:, end);
clear r_X;

r_dG0_f = r_S \ r_dG0_r;

figure(2);
plot(r_S * r_dG0_f, r_dG0_r, '.');