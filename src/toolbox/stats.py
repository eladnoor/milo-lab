import scipy.stats as st

def GetSlopeInterval(x, y, confidence=0.95):
    slope, _intercept, _r_value, _p_value, std_err = st.linregress(x,y)
    df = len(x) - 2 # unexplained degrees of freedom
    
    t = st.t.ppf((confidence + 1.0)/2.0, df)
    return slope - t*std_err, slope + t*std_err

def test():
    x = [0,12.0,29.5,43.0,53.0,62.5,75.5,85.0,93.0]
    y = [8.98,8.14,6.67,6.08,5.90,5.83,4.68,4.20,3.72]
    print GetSlopeInterval(x, y)
    
if __name__ == "__main__":
    test()