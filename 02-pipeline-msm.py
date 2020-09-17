#!/usr/bin/env python
# coding: utf-8

# ## <center>MSM-building Calculations</center>

# ##### MUST SHARE KERNEL WITH NOTEBOOK `00-load-data.ipynb`

# #### TICA and VAMP
# 
# ------------------------

# In[18]:

exec(open("00-load-data.py").read())


da = Path("analyses")
if not da.is_dir():
    da.mkdir()

#print(colorama.Back.CYAN+"Starting MSM Pipeline with TICA")
#print("Creating an array of models")
#
## uses last for `nm`,`feat` from kernel connection!!!
#n_ticas = len(all_models[feat][nm]["tica"])
#
#for n,(feat,nm,dataset) in enumerate(aswa_tools.iter_models(all_models)):
#    dds = da / nm
#
#    dtica = dds / "tica"
#    ftica = dtica / "results.pyemma"
#    
#    dvamp = dds / "vamp"
#    fvamp = dvamp / "results.pyemma"
#    
#    if not dds.is_dir():
#        dds.mkdir()
#    if not dtica.is_dir():
#        dtica.mkdir()
#    if not dvamp.is_dir():
#        dvamp.mkdir()
#
#    for m,setup in enumerate(dataset["tica"]):
#        
#        thisone = _model_name(feat, setup["kwargs"])
#        
#        _ftica = ftica.with_suffix(".%s"%thisone)
#        
#        if _ftica.is_file():
#            setup["result"] = pyemma.load(_ftica)
#            
#        else:
#            try:
#                setup["result"] = coor.tica(setup["input"], **setup["kwargs"])
#                setup["result"].save(_ftica, save_streaming_chain=True)
#
#            except:
#                print(colorama.Back.YELLOW+"Calculation failed!")
#                print("Saving the error in place of model: SOON! FIXME")
#                print(colorama.Back.RESET+colorama.Fore.LIGHTRED_EX+                    traceback.format_exc())
#
##    for setup in dataset["vamp"]:
# #       _fvamp = fvamp.with_suffix(".%s"%_model_name(feat, setup["kwargs"]))
#  #      
#   #     if _fvamp.is_file():
#    #        print("Found Results File: {}".format(_fvamp))
#     #       setup["result"] = pyemma.load(_fvamp)
#      #      
#       # else:
#        #    print("Calculating and Saving to File: {}".format(_fvamp))
#         #   setup["result"] = coor.vamp(setup["input"], **setup["kwargs"])
#          #  setup["result"].save(_fvamp, save_streaming_chain=True)
#
#
## #### k-means and regular-space clusterings
## ------------------------
#
## In[ ]:


# uses last for `nm`,`feat` from kernel connection!!!
n_kmeans    = len(all_models[feat][nm]["cluster_kmeans"])
n_regspaces = len(all_models[feat][nm]["cluster_regspace"])

for n,(feat,nm,dataset) in enumerate(aswa_tools.iter_models(all_models)):

    dds       = da / nm
    dkmeans   = dds / "kmeans" / ("tica_lag_%d"%chosen_lag)
    dregspace = dds / "regspace" / ("tica_lag_%d"%chosen_lag)
    fkmeans   = dkmeans / "results.pyemma"
    fregspace = dregspace / "results.pyemma"

    if not dkmeans.is_dir():
        dkmeans.mkdir(parents=True)
    if not dregspace.is_dir():
        dregspace.mkdir(parents=True)
    
    for m,setup in enumerate(dataset["cluster_kmeans"]):
        
        thisone = _model_name(feat, setup["kwargs"])
        _fkmeans = fkmeans.with_suffix(".%s"%thisone)
        
        if _fkmeans.is_file():
            print("Found Results File: {}".format(_fkmeans))
            setup["result"] = pyemma.load(_fkmeans)
            
        else:
            print("Calculating and Saving to File: {}".format(_fkmeans))
            
            if not setup["input"]:
                setup["input"] = aswa_tools.get_matching_input(
                    dataset["tica"], key="lag", val=chosen_lag)
                
                # SUM TING WONG with TICA object
                setup["input"]._default_chunksize = 10000
            
            tica_inp = np.array(np.concatenate(
                setup["input"].get_output(
                    list(range(n_tica_dim)))
                ), dtype=np.float64)
            
            setup["result"] = coor.cluster_kmeans(tica_inp, **setup["kwargs"])
            setup["result"].save(_fkmeans, save_streaming_chain=True)

    for m,setup in enumerate(dataset["cluster_regspace"]):

        thisone = _model_name(feat, setup["kwargs"])
        _fregspace = fregspace.with_suffix(".%s"%thisone)

       # if not setup["input"]:
       #     setup["input"] = aswa_tools.get_matching_input(
       #         dataset["tica"], key="lag", val=chosen_lag)
       #     # SUM TING WONG with TICA object
       #     setup["input"]._default_chunksize = 10000
       #         
       # tica_inp = np.array(np.concatenate(
       #     setup["input"].get_output(
       #         list(range(n_tica_dim)))
       #     ), dtype=np.float64)
       #     
       # mindist = np.power(np.product([
       #     mx-mn for mn,mx in zip(
       #         np.min(tica_inp[:,:n_tica_dim], axis=0),
       #         np.max(tica_inp[:,:n_tica_dim], axis=0))
       # ]), 1./n_tica_dim) / 25.
       # 
       # setup["actual-dmin"] = setup["kwargs"]["dmin"] * mindist
            
        if _fregspace.is_file():
            print("Found Results File: {}".format(_fregspace))
            setup["result"] = pyemma.load(_fregspace)
            
        else:
            print("Calculating and Saving to File: {}".format(_fregspace))
            setup["result"] = aswa_tools.refined_regspace(tica_inp, mindist, **setup["kwargs"])
            setup["result"].save(_fregspace, save_streaming_chain=True)


#for n,(feat,nm,dataset) in enumerate(aswa_tools.iter_models(all_models)):
#    
#    dataset["its"] = dlits = dict()
#
#    for m,(clust_method, par_key) in enumerate(
#        [("regspace", "dmin"), ("kmeans", "k")]):    
#        dlits[clust_method] = lits = list()
#        
#        dds = da / nm
#        dmsm = dds / "msm" / ("tica_lag_%d"%chosen_lag)
#        fits = dmsm / "its-results.pyemma"
#        
#        if not dmsm.is_dir():
#            dmsm.mkdir(parents=True)
#        
#        for i,setup in enumerate(
#            dataset["cluster_%s"%clust_method]
#        ):
#            thisone = _model_name(feat, setup["kwargs"])
#                    
#            print(colorama.Back.CYAN+                "Now creating MSMs using clusterings"+                "as the discretization to state space")
#
#            print(colorama.Back.CYAN+                "Calculated %d models so far"%(m))
#            print((" - on to model %s"%thisone)+                colorama.Back.RESET)
#            print(colorama.Fore.LIGHTGREEN_EX+                " - number %d / %d (wrong 2nd number)"%(
#                n_kmeans*n+m,(n_kmeans+n_regspaces)*len(all_models)))
#
#            _fits = fits.with_suffix(
#                ".%s-%s"%(clust_method, thisone))
#            
#            dits           = dict()
#            dits["par"]    = (par_key, setup["kwargs"][par_key])
#            dits["input"]  = setup["result"].dtrajs
#            dits["kwargs"] = dict(
#                lags=tica_lags,
#                reversible=False,
#                nits=3,
#                errors="bayes",
#            )
#            
#            if not _fits.is_file():
#                print("Calculating and Saving to File: {}".format(_fits))
#                dits["result"] = pyemma.msm.its(
#                    dits["input"],
#                    **dits["kwargs"]
#                )
#                dits["result"].save(_fits)
#            
#            else:
#                print("Found Results File: {}".format(_fits))
#                dits["result"] = pyemma.load(_fits)
#            
#            # Reverse so both with increasing N clusters
#            if "regspace" in clust_method:
#                lits.insert(0, dits)
#            else:
#                lits.append(dits)



reversible=False
for feat,nm,dataset in aswa_tools.iter_models(all_models):

    dataset["hits"] = dlits = dict()

    for clust_method, par_key in [("regspace", "dmin"), ("kmeans", "k")]:
    #for clust_method, par_key in [("kmeans", "k")]:
    #for clust_method, par_key in [("regspace", "dmin")]:
    
        dlits[clust_method] = lits = list()

        dds = da / nm
        dmsm = dds / "hmm" / ("tica_lag_%d"%chosen_lag)
        fits = dmsm / "its-results.pyemma"

        if not dmsm.is_dir():
            dmsm.mkdir(parents=True)

        for n in n_macrostates:
            for i,setup in enumerate(dataset["cluster_%s"%clust_method]):

                _fits = fits.with_suffix(
                    ".%s-%s_%d-%s_%s-%s"%(
                        clust_method, "nstates", n,
                        "reversible", str(reversible),
                        _model_name(feat, setup["kwargs"])))

                dits           = dict()
                dits["input"]  = setup["result"].dtrajs
                dits["par"]    = ((par_key, setup["kwargs"][par_key]),
                                  ("nstates", n))

                dits["kwargs"] = dict(
                    lags=tica_lags,
                    nstates=n,
                    reversible=reversible,
                    #nits=3,
                    errors="bayes",
                )

                if not _fits.is_file():
                    print("Calculating and Saving to File: {}".format(_fits))
                    pyemma.msm.timescales_hmsm(
                    #dits["result"] = pyemma.msm.timescales_hmsm(
                        dits["input"],
                        **dits["kwargs"]
                    #)
                    #dits["result"].save(_fits)
                    ).save(_fits)

                else:
                    print("Found Results File: {}".format(_fits))
                    #dits["result"] = pyemma.load(_fits)

        #        # Reverse so both with increasing N clusters
         #       if "regspace" in clust_method:
          #          lits.insert(0, dits)
           #     else:
            #        lits.append(dits)


# In[ ]:




