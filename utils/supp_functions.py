####################################################################################################
#                                        supp_functions.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: Aaron T. Gudmundson                                                                     #
#                                                                                                  #
# Created: 03/02/23                                                                                #
#                                                                                                  #
# Purpose: Python functions to perform meta-analysis operations. Used in the following work:       #
#                                                                                                  #
#           Gudmundson AT, Koo A, Virovka A, et al.                                                #
#           Meta-analysis and open-source database for in vivo brain Magnetic Resonance            #
#           spectroscopy in health and disease.                                                    #
#           Analytical Biochemistry. 2023;676:115227. doi:10.1016/j.ab.2023.115227                 # 
#                                                                                                  #
####################################################################################################
from scipy.stats import norm, ttest_ind, t
import numpy as np

def convert_median(n, median, q1, q3):
    '''
    Description:
      - Convert Median and Quartile Values to Mean and Standard Deviation
          Normal distributions uses Method described in Wan, 2014   (10.1186/1471-2288-14-135)
          Skewed distributions uses Method described in Greco, 2015 (10.13105/wjma.v3.i5.215 )
     
          *note Normal distributed data is defined as .5 > q1/q2 > 1.5.
    
    Inputs:
      - n      - Sample Size
      - median - Median Value
      - q1     - 1st Quartile (25%)
      - q3     - 3rd Quartile (75%)
      
    Outputs:
      - Results  -  Dictionary;
        - u      - Estimated mean
        - Std    - Estimated Standard Deviation
        - Ratio  - Ratio of Q1/Q3
        - Dist   - Estimated Distrubtion
        - Method - Method used for Estimation
    '''

    median     = float(median)
    q1         = float(q1)
    q3         = float(q3)
    
    Cochrane   = (q3 - q1)/1.35
    Wan_2014   = (q3 - q1) / (2 * norm.ppf((.75 * n - .125) / (n + .25)))
    Greco_2015 = np.array([(q3 - median)/.6745, (median - q1)/.6745])
    Greco_2015 = np.max(Greco_2015)

    ratio      = ((q3 - median) / (median - q1))
    
    print('Ratio   |  Med   |  Cochrane  |  Wan 2014  |  Greco 2015')

    if ratio < .5 or ratio > 1.5:
        print('{:.3f}    {:6.3f}     {:.3f}      {:.3f}     {:.3f}'.format(ratio, median, Cochrane, Wan_2014, Greco_2015))
        print('{:.3f}    {:6.3f}     {:.3f}      Greco 2015\n'.format(ratio, median, Greco_2015))

        u      = median
        std    = Greco_2015
        dist   = 'Skewed'
        meth   = 'Greco_2015'

    else:
        print('{:.3f}    {:6.3f}     {:.3f}      {:.3f}     {:.3f}'.format( ratio, median, Cochrane, Wan_2014, np.max(Greco_2015)))
        print('{:.3f}    {:6.3f}     {:.3f}      Wan 2014\n'.format(ratio, median, Wan_2014))

        u      = median
        std    = Wan_2014
        dist   = 'Normal'
        meth   = 'Wan_2014'

    results = {'u'     : u    ,
               'Std'   : std  ,
               'Ratio' : ratio,
               'Dist'  : dist ,
               'Method': meth }

    return results

def ratio_Friedrich(n1, n2, u1, u2, std1, std2):
    '''
    Description:
      - Determine the ratio from the mean and standard deviation of 2 groups using:
          - Friedrich, 2008 (10.1186/1471-2288-8-32)
          - Friedrich, 2011 (10.1016/j.jclinepi.2010.09.016)
         
    Inputs:
      - n1    - Group 1 Sample Size
      - n2    - Group 2 Sample Size
      - u1    - Group 1 Mean
      - u2    - Group 2 Mean
      - std1  - Group 1 Standard Deviation
      - std2  - Group 2 Standard Deviation
      
    Outputs:
      - Results  - Dictionary;
        - Ratio  - Estimated mean
        - Var    - Estimated Variance
        - SE     - Estimated Standard Error
        - Low    - Estimated Confidence Interval Low  (95%)
        - Hgh    - Estimated Confidence Interval High (95%)
    '''

    ratio     = (u1/u2)
    ln_ratio  = np.log(ratio)
    
    se0       = ((1/n1) * ((std1/u1)**2) )
    se1       = ((1/n2) * ((std2/u2)**2) )
    Variance  = (se0 + se1)
    SE        = np.sqrt( se0 + se1 )
    
    CI_low    = np.exp(ln_ratio - 1.96 * SE)
    CI_high   = np.exp(ln_ratio + 1.96 * SE)
    
    results   = {'Ratio': ratio,
                 'Var'  : Variance,
                 'SE'   : SE,
                 'Low'  : CI_low,
                 'Hgh'  : CI_high}

    return results

def effect_size(n1, n2, u1, u2, std1, std2):
    '''
    Description:
      - Calculate the Effect Size between 2 Groups using:
        - Cohen's D - Cohen, 1988 (ISBN-10: 0805802835; ISBN-13 978-0805802832)
        - Hedges' G - Hedge's, 1981 (10.2307/1164588)
         
    Inputs:
      - n1      - Group 1 Sample Size
      - n2      - Group 2 Sample Size
      - u1      - Group 1 Mean
      - u2      - Group 2 Mean
      - std1    - Group 1 Standard Deviation
      - std2    - Group 2 Standard Deviation
      
    Outputs:
      - Results - Dictionary;
        - Std_Pooled - Pooled Standard Deviation
        - CohensD    - Cohen's D
        - HedgesG    - Hedge's G (10.2307/1164588)
    '''
    
    u        = u1 - u2
    std1     = (n1 - 1) * (std1**2)
    std2     = (n2 - 1) * (std2**2)

    sd_pool  = (std1 + std2) / (n1 + n2 - 2) 
    sd_pool  = np.sqrt(sd_pool)
        
    d        = u/sd_pool    
    g        = d * (1 - (3 / ((4 * (n1 + n2)) - 9)))

    results  = {'Std_Pooled': sd_pool,
                'CohensD'   : d      ,
                'HedgesG'   : g      }

    return results
    

def fixed_effect(df, Tcol, Wcol, x0=0.0):
    '''
    Description:
      - Calculate the Combined Random Effect using the methods described:
          - Borenstein, 2009 (10.1002/9780470743386)
          - Borenstein, 2010 (10.1002/jrsm.12)
         
    Inputs:
      - df      - Pandas DataFrame; All Studies and respective Mean and Weights
      - Tcol    - String; Mean Column
      - Tcol    - String; Weight Column
      - x0      - Float;  Specific Value to Test Null-Hypothesis
      
    Outputs:
      - Results - Dictionary:
          - M          - Weighted Mean Value
          - Var        - Variance
          - Weight     - Sum of Weights
          - SE         - Standard Error
          - Z          - Z (Z-distrubution)
          - P-1Tail    - P-value (1-way)
          - P-2Tail    - P-value (2-way)
          - Low        - Estimated Confidence Interval Low  (95%)
          - Hgh        - Estimated Confidence Interval High (95%)
    '''    

    Tcol    = df[Tcol].values                                                       # 
    Wcol    = df[Wcol].values
    
    M       = Tcol * Wcol
    M_u     = np.sum(M)/np.sum(Wcol)
    V       = 1/np.sum(Wcol)
    SE      = np.sqrt(V)
    Low     = M_u - (1.96 * SE)
    Hgh     = M_u + (1.96 * SE)
    
    Z       = (M_u - x0)/SE                                                         # Z; Z-distribution
    p1way   = 1 - norm.cdf(np.abs(Z))
    p2way   = 2 * p1way
    
    results = {'M'      : M_u,
               'Var'    : V  ,
               'SE'     : SE ,
               'Z'      : Z  ,
               'P-1Tail': p1way,
               'P-2Tail': p2way,
               'Weight' : np.sum(Wcol),
               'Low'    : Low,
               'Hgh'    : Hgh}

    return results
    
def random_effect(df, Tcol, Wcol, x0=0.0):
    '''
    Description:
      - Calculate the Combined Random Effect using the methods described:
          - Borenstein, 2009 (10.1002/9780470743386)
          - Borenstein, 2010 (10.1002/jrsm.12)

    Inputs:
      - df      - Pandas DataFrame; All Studies and respective Mean and Weights
      - Tcol    - String; Mean Column
      - Tcol    - String; Weight Column
      - x0      - Float;  Specific Value to Test Null-Hypothesis
     
    Outputs:
      - Results - Dictionary; With Random Effects Model
          - M          - Weighted Mean Value
          - Var        - Variance
          - Weight     - Sum of Weights
          - SE         - Standard Error
          - Q          - Weighted Sum of Squares
          - Z          - Z (Z-distrubution)
          - P-1Tail    - P-value (1-way)
          - P-2Tail    - P-value (2-way)
          - Low        - Estimated Confidence Interval Low  (95%)
          - Hgh        - Estimated Confidence Interval High (95%)
    '''    
    
    Tcol    = df[Tcol].values                                                       # Mean Values
    Wcol    = df[Wcol].values                                                       # Weight Values
  
    ## Between Study Variance -  DerSimonian and Laird Method (10.1016/0197-2456(86)90046-2)
    #  T2 = (Q - df) / C

    ## Q = Weighted Sum of Squares
    Q1      = (Tcol ** 2) * Wcol                                                    # Weighted Sum of Squares (Left)
    Q1      = np.sum(Q1)                                                            # Weighted Sum of Squares (Left)
    
    Q2      = np.sum(Tcol * Wcol)                                                   # Weighted Sum of Squares (Top)
    Q2      = Q2**2                                                                 # Weighted Sum of Squares (Top)
    
    Q3      = np.sum(Wcol)                                                          # Weighted Sum of Squares (Bottom)
    
    Q       = Q1 - (Q2/Q3)                                                          # Weighted Sum of Squares (Full)
 
    ## df = Degrees of Freedom
    dfrdm   = df.shape[0] - 1                                                       # Degrees of Freedom

    I       = (Q - df.shape[0])/Q                                                   # Heterogeneity
    I      *= 100                                                                   # Heterogeneity
    
    ## C = Constant 
    #   "Since Q - df is on standardized scale... puts index back into the same metric that had been used to report the within study variance."
    #      - From Borenstein, 2010 (10.1002/jrsm.12)
    C1      = np.sum(Wcol)                                                          # Between Study Variance 
    C3      = np.sum(Wcol)                                                          # Between Study Variance 
    
    C2      = Wcol ** 2                                                             # Between Study Variance
    C2      = np.sum(C2)                                                            # Between Study Variance
    C       = C1 - (C2/C3)                                                          # Between Study Variance
        
    ## T2
    Tau_2   = (Q - dfrdm)/C                                                         # Between Study Variance - Tau Squared

    if Tau_2 < 0:                                                                   # If Tau Squared is Negative         
        Tau_2 = 0                                                                   # Set to 0       

    Wcol    = 1/Wcol                                                                # Weight   --> Variance
    Wcol   += Tau_2                                                                 # Within Study Variance                                                            
    Wcol    = 1/Wcol                                                                # Variance --> Weights
    
    M       = Tcol * Wcol                                                           # Summary Effect Size
    M_u     = np.sum(M)/np.sum(Wcol)                                                # Summary Effect Size
        
    V       = 1/np.sum(Wcol)                                                        # Variance
    SE      = np.sqrt(V)                                                            # Standard Error
    Low     = M_u - (1.96 * SE)                                                     # Confidence Interval - Low
    Hgh     = M_u + (1.96 * SE)                                                     # Confidence Interval - High
    
    Z       = (M_u - x0)/SE                                                         # Z; Z-distribution
    p1way   = 1 - norm.cdf(np.abs(Z))                                               # 1-Tailed p-Value
    p2way   = 2 * p1way                                                             # 2-Tailed p-Value
    
    results = {'M'      : M_u,
               'Var'    : V  ,
               'SE'     : SE ,
               'Z'      : Z  ,
               'Q'      : Q  ,
               'I2'     : I  ,
               'P-1Tail': p1way,
               'P-2Tail': p2way,
               'Weight' : np.sum(Wcol),
               'Low'    : Low,
               'Hgh'    : Hgh}

    return results

def confidence(d, se):
    '''
    Description:
      - Calculate the Confidence Intervals using Effect Size:

    Inputs:
      - n1      - Sample Size Group 1
      - n2      - Sample Size Group 2
      - d       - Pandas DataFrame; All Studies and respective Mean and Weights
    
    Outputs:
      - results - dictionary;
          - Low - (Confidence Interval (Low )
          - Hgh - (Confidence Interval (High)

    '''

    dlow     = d - (1.96 * se)
    dhgh     = d + (1.96 * se)
    results  = {'Low': dlow,
                'Hgh': dhgh}
    return results

def conf_from_stats(u, std, n):
    '''
    Description:
      - Calculate the Confidence Intervals from Summary Statistics

    Inputs:
      - u       - Mean 
      - std     - Standard Deviation
      - n       - Sample Size 
    
    Outputs:
      - results - dictionary;
          - Low - (Confidence Interval (Low )
          - Hgh - (Confidence Interval (High)
    '''

    if np.isnan(u) == True:
        return np.nan, np.nan
    
    if n < 30:
        crt = t.ppf( (1-(.05/2)), (n-1))
    else:
        crt = norm.ppf(1-.05/2, (n-1))
    
    CI     = crt * (std/np.sqrt(n))
    CI_Low = u - CI
    CI_Hgh = u + CI
    return CI_Low, CI_Hgh
