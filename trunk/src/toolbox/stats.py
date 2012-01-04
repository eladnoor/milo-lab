import numpy as np
import scipy.stats as st

def MeanWithConfidenceInterval(Y, confidence=0.95):
    """
    Use the fact that (mean(Y) - mu) / (std(Y)/sqrt(n))
    is a Student T distribution with n-1 degrees of freedom
    
    Returns:
        2 tuple (mean, symmetric confidence interval size).
    """
    n = len(Y)
    Y_bar = st.nanmean(Y)

    # According to the Student T-test distribution for n-1 degrees of freedom
    # find the position where the CDF is 0.975 (assuming we want a confidence
    # of 0.95). The lower part of the tail will account for the other 0.025
    # chance.
    t = st.t.ppf((confidence + 1.0)/2.0, n-1)
    SD = st.nanstd(Y, bias=False) # use the unbiased estimator: sqrt(y^2 / (n-1))
    SE = SD / np.sqrt(len(Y))
    return Y_bar, t*SE


def MeanWithInterval(Y, confidence=0.95):
    """
        Use the fact that (mean(Y) - mu) / (std(Y)/sqrt(n))
        is a Student T distribution with n-1 degrees of freedom
    """
    Y_bar, err = MeanWithConfidenceInterval(Y, confidence)
    
    return Y_bar - err, Y_bar + err




def GetSlopeInterval(x, y, confidence=0.95):
    slope, _intercept, _r_value, _p_value, std_err = st.linregress(x,y)
    df = len(x) - 2 # unexplained degrees of freedom
    
    t = st.t.ppf((confidence + 1.0)/2.0, df)
    return slope - t*std_err, slope + t*std_err

def test():
    x = [0,12.0,29.5,43.0,53.0,62.5,75.5,85.0,93.0]
    y = [8.98,8.14,6.67,6.08,5.90,5.83,4.68,4.20,3.72]
    
    print MeanWithInterval(x)
    print MeanWithInterval(y)
    print GetSlopeInterval(x, y)
    
if __name__ == "__main__":
    test()