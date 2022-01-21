#%%
import numpy as np
import random

#%%
s_num = 25
t_num = 25
pt_lb = 1
pt_ub = 25
pre_lb = 1
pre_ub = 25

# %%
def generate_suc_by_pre(pre,s_num,t_num):
    suc=[[] for i in range(s_num)]
    for i in range(t_num):
        for s in pre[i]:
            suc[s].append(i)
    return suc
def generate_data(s_num,t_num,pt_lb,pt_ub,pre_lb,pre_ub):#processing time上下界，前置job上下界
    p_time_setup=[random.randrange(pt_lb,pt_ub) for i in range(s_num)]
    p_time_task=[random.randrange(pt_lb,pt_ub) for i in range(t_num)]
    pre=[random.sample(range(s_num), random.randrange(pre_lb,pre_ub)) for i in range(t_num)]
    suc=generate_suc_by_pre(pre,s_num,t_num)
    for s in range(s_num):
        if suc[s]==[]:
            t=random.randrange(t_num)
            suc[s].append(t)
            pre[t].append(s)
    return pre,suc,p_time_setup,p_time_task
def tran_to_graph(suc):
    g=np.zeros((s_num,t_num))
    for i in range(s_num):
        for j in suc[i]:
            g[i][j]=1
    # print(g)
    for i in range(len(g)):
        for j in range(len(g[0])):
            print(int(g[i][j]),end=' ')
        print()
#%%
pre,suc,p_time_setup,p_time_task=generate_data(s_num,t_num,pt_lb,pt_ub,pre_lb,pre_ub)
s_num,t_num=len(p_time_setup),len(p_time_task)

# %%
print(f"{s_num} {t_num}")
for i in range(len(p_time_setup)):
    print(f"{p_time_setup[i]}", end=" ")
print()
for i in range(len(p_time_task)):
    print(f"{p_time_task[i]}", end=" ")
print()
print("-1")
tran_to_graph(suc)

#%%
# %%
