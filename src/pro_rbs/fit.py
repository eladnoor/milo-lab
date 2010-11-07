import csv, os, sys
from pylab import *
from scipy.optimize import leastsq
from toolbox.html_writer import HtmlWriter
from toolbox import util

"""
    Read the data from the CSV files containing measurements of the activity of a gene
    as a function of the RBS (row) and promoter (columns)
"""

def read_data(fname):
    facs_csv = csv.reader(open(fname, 'r'))
    promoters = facs_csv.next()[1:] # the first column is the RBS
    rbs = []
    data = []
    for row in facs_csv:
        rbs.append(row[0])
        data.append([float(x) for x in row[1:]])

    data = matrix(data)
    return (rbs, promoters, data)

def fit_err(param, data, threshold):
    err = []
    [Nr, Np] = data.shape
    R = param[:Nr]
    P = param[Nr:(Nr+Np)]
    floor = data.min()
    for i in range(Nr):
        for j in range(Np):
            if (data[i, j] >= threshold):
                err.append(R[i] + P[j] - log10(data[i, j]))
            else:
                err.append(0)
    return array(err)

def fit_err2(param, data):
    err = []
    [Nr, Np] = data.shape
    R = param[:Nr]
    P = param[Nr:(Nr+Np)]
    C = param[Nr+Np:]
    for i in range(Nr):
        for j in range(Np):
            err.append(log(data[i, j]) - log(C[0] + 10**(R[i] + P[j])))
    return array(err)

def fit_params(csv_fname, html, R_calc, threshold=None):
    (rbs, promoters, data) = read_data('../data/pro_rbs/' + csv_fname + ".csv")
    [Nr, Np] = data.shape
    if (threshold == None): # fit this function: C + 10^(R_i + P_j) -- in logscale
        x0 = zeros(Nr+Np+1)
        (param, ier) = leastsq( (lambda x:fit_err2(x, data)), x0)
        C = param[Nr+Np:]
    else: # fit this function: 10^(R_i + P_j) -- in logscale, ignore examples which are below the threshold
        x0 = zeros(Nr+Np)
        (param, ier) = leastsq( (lambda x:fit_err(x, data, threshold)), x0)
        C = [0]

    R = param[:Nr]
    P = param[Nr:(Nr+Np)]
    b = min(R) + min(P)
    R = R - R.min()
    P = P - P.min()
    est = []
    mes = []
    for i in range(Nr):
        for j in range(Np):
            mes.append(data[i, j])
            est.append(C[0] + 10**(b + R[i] + P[j]))

    est = array(est)
    mes = array(mes)
    
    fig1 = figure()
    hold(True)
    if (threshold == None):
        plot([mes.min(), mes.max()], [mes.min(), mes.max()], 'k--')
        plot([C[0], C[0]], [mes.min(), mes.max()], 'b:')
        plot(est, mes, '.g')
        r2 = corrcoef(est, mes)[0,1]**2
    else:
        plot([threshold, mes.max()], [threshold, mes.max()], 'k--')
        plot([threshold, threshold], [mes.min(), mes.max()], 'b:')
        plot(est[find(mes >= threshold)], mes[find(mes >= threshold)], '.g')
        plot(threshold*ones(len(find(mes < threshold))), mes[find(mes < threshold)], '.b')
        r2 = corrcoef(est[find(mes >= threshold)], mes[find(mes >= threshold)])[0,1]**2
    title(csv_fname + ', r^2 = %.2f' % r2)
    xscale('log')
    yscale('log')
    xlabel('Estimation')
    ylabel('Measurement')

    html.write('Global offset = %.1f<br>\n' % b)
    html.write('<table border="1">\n')
    html.write('<tr><td>RBS\promoter</td>')
    for j in range(Np):
        html.write('<td>%s (%.1f)</td>' % (promoters[j], P[j]))
    html.write('</tr>\n')

    for i in range(Nr):
        html.write('<tr>')
        html.write('<td>%s (%.1f)</td>' % (rbs[i], R[i]))
        for j in range(Np):
            mes = log10(data[i, j])
            est = log10(C[0] + 10**(b + R[i] + P[j]))
            if (threshold != None and mes < log10(threshold)):
                html.write('<td><font color="blue">%.1f</font></td>' % mes)
            elif (abs(mes-est) > 0.3):
                html.write('<td><font color="red">%.1f (%.1f)</font></td>' % (mes, est))
            else:
                html.write('<td>%.1f (%.1f)</td>' % (mes, est))
        html.write('</tr>\n')
    html.write('</table>\n')
    html.write('Measured (Estimated) - in logarithmic scale (base 10)</br>\n')
    html.write('red marks estimations that are off by more than 0.3</br>\n')
    if (threshold != None):
        html.write('blue marks estimations that are below the threshold (%.1f)</br>\n' % log10(threshold)  )
    
    savefig('../res/pro_rbs/%s_correlation.pdf' % csv_fname, format='pdf')
    html.embed_matplotlib_figure(fig1, width=800, height=600)

    if (R_calc != None):
        figure()
        plot(R_calc, 10**R, '.')
        xscale('log')
        yscale('log')
        xlabel('RBS strength (Voigt)')
        ylabel('RBS strength (current model)')
        
        savefig('../res/pro_rbs/%s_rbs_compare.svg' % 'data_flu-OD_multi', format='svg')
        html.embed_svg('%s_rbs_compare.svg' % 'data_flu-OD_multi', width=800, height=600)

    return (b, R, P)

def read_rbs_calc_results(fname):
    csv_reader = csv.reader(open(fname, 'r'))
    csv_reader.next()
    R_calc = []
    for row in csv_reader:
        (Name,Start_Position,Expression_Level,Kinetic_Score,Sequence) = row
        R_calc.append(float(Expression_Level))
    return array(R_calc)

################################################################################
################################################################################
################################################################################
################################################################################

#(rbs, promoters, data_facs) = read_data('data_FACS.csv')
#(rbs, promoters, data_single) = read_data('data_flu-OD_single.csv')
#(rbs, promoters, data_multi) = read_data('data_flu-OD_multi.csv')

util._mkdir('../res/pro_rbs')

rbs_score_fname = '../res/pro_rbs/rbs_2010-08-18_17-50-19_133.csv'
if (os.path.exists(rbs_score_fname)):
    R_calc = read_rbs_calc_results(rbs_score_fname)
else:
    sys.stderr.write("The RBS calculator score file could not be found, you must " + \
                    "generate it using the 'rbs-calc' website and put it here: %s" % rbs_score_fname)
    R_calc = None
    
html = HtmlWriter('../res/pro_rbs/fit.html')
(b, R, P) = fit_params('data_flu-OD_multi', html, R_calc)
