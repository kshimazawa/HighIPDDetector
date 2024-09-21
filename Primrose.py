import pysam
import numpy as np
import sys
import matplotlib.pyplot as plt
import lightgbm as lgb
import pandas as pd
import shap
import seaborn as sns
def primrose(RESULT_PATH):
    # kinetics = []
    # context_list = []
    # ml_list = []
    # ml_cnt = 0
    # bamfile = pysam.AlignmentFile(ALIGNED_FILE,'rb',check_sq=False)
    # for read in bamfile.fetch(until_eof=True):
    #     mean_ipd_f = np.array(read.get_tag('fi'),dtype=int)
    #     mean_ipd_r = np.flip(np.array(read.get_tag('ri'),dtype=int))
    #     mean_pw_f = np.array(read.get_tag('fp'),dtype=int)
    #     mean_pw_r = np.flip(np.array(read.get_tag('rp'),dtype=int))
    #     if (len(mean_ipd_f)==0) | (len(mean_ipd_r)==0):
    #         print("no values for either forward/reverse strand")
    #         continue
    #     try:
    #         ml = np.array(read.get_tag('ML'),dtype=int)
    #     except:
    #         print("no ml tag")
    #         continue
    #     mm = np.array(read.get_tag('MM').strip(";").split(",")[1:],dtype=int)
    #     seq = read.get_forward_sequence()
    #     mm_iter = 0
    #     context = []
    #     c_iter = 0
    #     for i in range(len(seq)):
    #         if mm_iter == len(mm):
    #             break
    #         if (seq[i] != "C"):
    #             continue
    #         if c_iter == np.sum(mm[:mm_iter+1])+mm_iter:
    #             # if (seq[i-1:i+4] == "CCAGG") | (seq[i-1:i+4] == "CCTGG"):
    #             #     context.append(ml[mm_iter])
    #                 # print("appended")
    #             # if (i >= 8) & (i<=len(seq)-11) & (ml[mm_iter] >= 128):
    #             #     context.append(seq[i-8:i+11])
    #             mm_iter += 1
    #             f_start_iter = max(i-5,0)
    #             f_end_iter = min(i+11,len(seq))
    #             r_start_iter = max(i-8,0)
    #             r_end_iter = min(i + 8,len(seq))
    #             kinetics_row = np.zeros((4,16))
    #             if f_start_iter == 0:
    #                 kinetics_row[0][5-i:] = mean_ipd_f[f_start_iter:f_end_iter]
    #                 kinetics_row[2][5-i:] = mean_pw_f[f_start_iter:f_end_iter]
    #             elif f_end_iter == len(seq):
    #                 kinetics_row[0][:(f_end_iter-f_start_iter)] = mean_ipd_f[f_start_iter:f_end_iter]
    #                 kinetics_row[2][:(f_end_iter-f_start_iter)] = mean_pw_f[f_start_iter:f_end_iter]
    #             else:
    #                 kinetics_row[0] = mean_ipd_f[f_start_iter:f_end_iter]
    #                 kinetics_row[2] = mean_pw_f[f_start_iter:f_end_iter]
    #             if r_start_iter == 0:
    #                 kinetics_row[1][8-i:] = mean_ipd_r[r_start_iter:r_end_iter]
    #                 kinetics_row[3][8-i:] = mean_pw_r[r_start_iter:r_end_iter]
    #             elif r_end_iter == len(seq):
    #                 kinetics_row[1][:(r_end_iter-r_start_iter)] = mean_ipd_r[r_start_iter:r_end_iter]
    #                 kinetics_row[3][:(r_end_iter-r_start_iter)] = mean_pw_r[r_start_iter:r_end_iter]
    #             else:
    #                 kinetics_row[1] = mean_ipd_r[r_start_iter:r_end_iter]
    #                 kinetics_row[3] = mean_pw_r[r_start_iter:r_end_iter]
    #             kinetics.append(kinetics_row)
    #         c_iter += 1
    #     ml_list.append(ml)
    #     ml_cnt += len(ml)
    #     context_list.append(context)
    # bamfile.close()
    # kinetics = np.array(kinetics)
    # ml_flatten = np.array([v for ml in ml_list for v in ml])
    # # category =["forward_ipd","reverse_ipd","forward_pw","reverse_pw"]
    # df_X = pd.DataFrame(np.reshape(kinetics,(kinetics.shape[0],kinetics.shape[1]*kinetics.shape[2])),columns=np.arange(kinetics.shape[1]*kinetics.shape[2]))
    # df_Y = pd.DataFrame(ml_flatten,columns=['y'])
    # fold = 10
    # print("learning started")
    # label_fi = ["fIPD{}".format(i) for i in range(-5,11)]
    # label_fp = ["fPW{}".format(i) for i in range(-5,11)]
    # label_ri = ["rIPD{}".format(i) for i in range(-8,8)]
    # label_rp = ["rPW{}".format(i) for i in range(-8,8)]
    # label_total = label_fi+label_ri+label_fp+label_rp
    df_me = pd.read_csv(RESULT_PATH+"methylated_sub.csv")
    df_un = pd.read_csv(RESULT_PATH+"unmethylated_sub.csv")
    df = pd.concat([df_me,df_un]).sample(frac=1,ignore_index=True)
    df_X = df.loc[:,'fIPD-5':'rPW7']
    df_Y = df.loc[:,'y']
    fold = 10
    for it in range(1):
        unit = len(df_X)//fold
        if it == 0:
            df_X_train =df_X[unit * 1:unit*fold]
            df_Y_train =df_Y[unit * 1:unit*fold]
        elif it == fold - 1:
            df_X_train =df_X[0:unit * it]
            df_Y_train =df_Y[0:unit * it]
        else:
            df_X_train =pd.concat([df_X[0:unit * it],df_X[unit * (it+1):unit*fold]])
            df_Y_train = pd.concat([df_Y[0:unit * it],df_Y[unit * (it+1):unit*fold]])
        df_X_test = df_X[unit * it: unit*(it+1)]
        df_Y_test = df_Y[unit * it: unit*(it+1)]
        lgb_train = lgb.Dataset(df_X_train,df_Y_train)
        lgb_eval = lgb.Dataset(df_X_test,df_Y_test)
        params = {
            'seed':4,
            'objective':'regression',
            'metric':'rmse'
        }
        evals_result = {}
        lgbm = lgb.train(params,lgb_train,valid_sets=[lgb_train,lgb_eval],valid_names=['train','eval'],num_boost_round=200,callbacks=[lgb.record_evaluation(evals_result),lgb.log_evaluation(period=100)])
        plt.plot(evals_result['train']['rmse'], label='train rmse')
        plt.plot(evals_result['eval']['rmse'], label='eval rmse')
        plt.legend()
        plt.xlabel('rounds')
        plt.ylabel('rmse')
        plt.savefig(RESULT_PATH+f"learning_primrose{it}.png")
        plt.close()
        explainer = shap.TreeExplainer(model=lgbm)
        shap_values = explainer.shap_values(X=df_X_test)
        coefs = np.zeros(64)
        for i in range(64):
            coef, intercept = np.polyfit(df_X_test.iloc[:,i].values,shap_values[:,i],1)
            coefs[i] = coef
        coefs = coefs.reshape((4,16))
        df_coefs = pd.DataFrame(np.concatenate([coefs[:1,1:11],coefs[1:2,4:14]]),index=["Forward IPD","Reverse IPD"],columns=list("NNNNCGNNNN"))
        plt.figure(figsize=(10,5))
        sns.heatmap(df_coefs,cmap='Oranges',vmin=-1.,vmax=1.,annot=True)
        plt.title("Coefficients of Linear Regression between IPDs and Their SHAP Values")
        plt.savefig("IPD.png")
        plt.close()
        plt.figure(figsize=(10,5))
        df_coefs = pd.DataFrame(np.concatenate([coefs[2:3,1:11],coefs[3:,4:14]]),index=["Forward PW","Reverse PW"],columns=list("NNNNCGNNNN"))
        plt.title("Coefficients between PWs and Their SHAP Values")
        sns.heatmap(df_coefs,cmap='Oranges',vmin=-4.,vmax=4.,annot=True)
        plt.savefig("PW.png")
        #shap.summary_plot(shap_values,df_X_test,feature_names=df_X.columns,show=False)
        #plt.savefig(RESULT_PATH+f"shap{it}.png")
        # lgb.plot_importance(lgbm,figsize=(10,10),max_num_features=10,importance_type='gain')
        # plt.savefig(RESULT_PATH+f"importance_features{it}.png")
        # plt.close()
    # for kin in range(4):
    #     for pos in range(16):
    #         plt.figure(figsize=(10,10))
    #         prob = []
    #         kinetics_to_use = kinetics[:,kin,pos]
    #         for i in range(10):
    #             prob.append(kinetics_to_use[(ml_flatten>256//10*i) & (ml_flatten<256//10*(i+1))])
    #         plt.boxplot(prob)
    #         # plt.hist2d(kinetics[:,kin,pos],ml_flatten,bins=[np.linspace(0,256,257),np.linspace(0,256,257)],density=True,cmap=cm.jet)
    #         plt.savefig(RESULT_PATH+"{}_pos{}.png".format(category[kin],pos))
    #         plt.close()
    # with open(RESULT_PATH+"context_fp.fasta","w") as fw:
    #     motif_cnt = 0
    #     for e in context_list:
    #         for m in e:
    #             fw.write(f">{motif_cnt}\n")
    #             fw.write(f"{m}\n")
    #             motif_cnt += 1
    return

        
if __name__ == "__main__":
    primrose(sys.argv[1])
