import statsmodels.api as sm
predata1 = []
for i in range(1,228):
    nn1 = 1
    n3 = 999+i
    n4=n3+72
    data1 = [rain1[n3-1:n4],water_lv1[n3-1:n4]]
    subset_data = water_lv1[n3-1:n4]
    model = sm.tsa.AutoReg(subset_data, lags=1)
    ar1 = model.fit()
    predata0 = ar1.forecast(steps=1)
    predata1+= [predata0]
predata1 = pd.Series(data=predata1,index = [i for i in range(len(predata1))])



