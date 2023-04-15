library(plyr)
library(dplyr) #����select
library(psych) #���ϵ������
library(pacman) #ȱʧֵ����
library(mice)
library("VIM")
library(pheatmap)
library(car)
library(ggplot2)
library(corrplot)
library(RColorBrewer)

######��������

setwd("d:/fx")
ms = read.table(file = "ms.csv", header = T , 
                     sep = "," , fill = TRUE , encoding = "UTF-8")
colnames(ms)[1] = "����"
head(ms)
ms = ms[1:45,]

######�޸�BACS�������Ϊzscore
ms$bacs���Լ��� = scale(ms$bacs���Լ���)
ms$bacs�������� = scale(ms$bacs��������)
ms$bacs���Ҷ��� = scale(ms$bacs���Ҷ���)
ms$bacs�׶��� = scale(ms$bacs�׶���)
ms$bacs������� = scale(ms$bacs�������)
ms$bacs���������� = scale(ms$bacs����������)
ms$bacs�ܷ� = scale(ms$bacs�ܷ�)

######BACS�����Բ���PANSS��סԺ��Ṧ�ܡ�Ѫ��ָ��������
results1 <- matrix(0,nrow=48,ncol=7,
                   dimnames=list(rep(c("t","p"),24),(colnames(ms)[c(24:30)])))
#summary(lm(bacs�ܷ�~��ϸ������+bmi+�Ա�+����,ms))
for (i in colnames(ms)[c(24:30)]){
  s = 0
  for (j in colnames(ms)[c(18:23,34:37,31,41:53)]){
    s = s+1
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms))
    results1[(s-1)*2+1,i] = fit$coefficients[2,3]
    results1[(s-1)*2+2,i] = fit$coefficients[2,4]
  }
}
write.csv(results1,"msresult1.csv")


######סԺ��Ṧ�������Բ���panss��Ѫ��������
results2 <- matrix(0,nrow=16,ncol=2,
                   dimnames=list(colnames(ms)[c(18:23,34:37,48:53)],c("t","p")))
for (i in colnames(ms)[c(18:23,34:37,48:53)]){
  formula = paste0("סԺ��Ṧ��~",i,"+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
    results2[i,1] = fit$coefficients[2,3]
    results2[i,2] = fit$coefficients[2,4]
}
write.csv(results2,"msresult2.csv")

######CRP��BACS�����Բ���PANSS��סԺ��Ṧ�ܡ�Ѫ��ָ��������
results3 <- matrix(0,nrow=25,ncol=2,
                   dimnames=list(colnames(ms)[c(18:32,34:37,48:53)],
                                 c("t","p")))
for (i in colnames(ms)[c(18:32,34:37,48:53)]){
  formula = paste0("crp~",i,"+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
  results3[i,1] = fit$coefficients[2,3]
  results3[i,2] = fit$coefficients[2,4]
}
write.csv(results3,"msresult3.csv")

######�Ƿ��м�״�ټ�������WHOQOL��סԺ����������WHO-DAS��BACS��PANSS��Ӱ��
ms$��״�ټ��� = as.factor(ms$��״�ټ���)
result3_5 <- matrix(0,nrow=20,ncol=3,
                  dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],c("coef","tstat","pvalue")))

summary(lm(bacs�ܷ�~��״�ټ���+bmi+�Ա�+����,ms))
for (i in colnames(ms)[c(9:15,24:32,34:37)]){
  formula = paste0(i,"~��״�ټ���+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
  result3_5[i,"coef"] = fit$coefficients[2,1]
  result3_5[i,"tstat"] = fit$coefficients[2,3]
  result3_5[i,"pvalue"] = fit$coefficients[2,4]
}

write.csv(result3_5,"msresult3_5.csv")



######����״�ټ��ض���WHOQOL��סԺ����������WHO-DAS��BACS��PANSS��Ӱ��
result4 <- matrix(0,nrow=20,ncol=3,
                     dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                                   c("coef","tstat","pvalue")))
ans = NA
#summary(lm(bacs�ܷ�~�����״����+bmi+�Ա�+����,ms))
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms))
    result4[i,"coef"] = fit$coefficients[2,1]
    result4[i,"tstat"] = fit$coefficients[2,3]
    result4[i,"pvalue"]  = fit$coefficients[2,4]
  }
  ans = cbind(ans,result4)
}
write.csv(ans,"msresult4.csv")

#result4���ϵ����ͼ
vis4 = matrix(0,nrow=20,ncol=5,
              dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                            colnames(ms)[c(43:47)]))
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms))
    vis4[i,j] = fit$coefficients[2,4]
  }
}
# ����Ԫ��ĸ߶ȺͿ��ȶ�����Ϊ20
pheatmap(-log(vis4,base=10), border_color = "black", #�߿���Ϊ��ɫ
         display_numbers = F,         #��ͼ��������ʾ��Ӧ����ֵ
         number_color = "black",         #������ɫΪ��ɫ
         fontsize=10,                    #�����СΪ10
         cluster_rows=F,                 #�����о���
         cluster_cols = F,
         number_format = "%.3f",         #����һλС��
         fontface="italic")              #����������Ϊб��

######(���м�״�ټ���)����״�ټ��ض���WHOQOL��סԺ����������WHO-DAS��BACS��PANSS��Ӱ��
result5 <- matrix(0,nrow=20,ncol=3,
                  dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                                c("coef","tstat","pvalue")))
ans = NA
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms[which(ms$��״�ټ���==1),]))
    result5[i,"coef"] = fit$coefficients[2,1]
    result5[i,"tstat"] = fit$coefficients[2,3]
    result5[i,"pvalue"]  = fit$coefficients[2,4]
  }
  ans = cbind(ans,result5)
}
write.csv(ans,"msresult5.csv")

#result5���ϵ����ͼ
vis5 = matrix(0,nrow=20,ncol=5,
              dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                            colnames(ms)[c(43:47)]))
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms[which(ms$��״�ټ���==1),]))
    vis5[i,j] = fit$coefficients[2,4]
  }
}
# ����Ԫ��ĸ߶ȺͿ��ȶ�����Ϊ20
pheatmap(-log(vis5,base=10), border_color = "black", #�߿���Ϊ��ɫ
         display_numbers = F,         #��ͼ��������ʾ��Ӧ����ֵ
         number_color = "black",         #������ɫΪ��ɫ
         fontsize=10,                    #�����СΪ10
         cluster_rows=F,                 #�����о���
         cluster_cols = F,
         number_format = "%.3f",         #����һλС��
         fontface="italic")              #����������Ϊб��

######(δ���м�״�ټ���)����״�ټ��ض���WHOQOL��סԺ����������WHO-DAS��BACS��PANSS��Ӱ��
result6 <- matrix(0,nrow=20,ncol=3,
                  dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                                c("coef","tstat","pvalue")))
ans = NA
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms[which(ms$��״�ټ���==0),]))
    result6[i,"coef"] = fit$coefficients[2,1]
    result6[i,"tstat"] = fit$coefficients[2,3]
    result6[i,"pvalue"]  = fit$coefficients[2,4]
  }
  ans = cbind(ans,result6)
}
write.csv(ans,"msresult6.csv")
#result6���ϵ����ͼ
vis6 = matrix(0,nrow=20,ncol=5,
              dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                            colnames(ms)[c(43:47)]))
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms[which(ms$��״�ټ���==0),]))
    vis6[i,j] = fit$coefficients[2,4]
  }
}
# ����Ԫ��ĸ߶ȺͿ��ȶ�����Ϊ20
pheatmap(-log(vis6,base=10), border_color = "black", #�߿���Ϊ��ɫ
         display_numbers = F,         #��ͼ��������ʾ��Ӧ����ֵ
         number_color = "black",         #������ɫΪ��ɫ
         fontsize=10,                    #�����СΪ10
         cluster_rows=F,                 #�����о���
         cluster_cols = F,
         number_format = "%.3f",         #����һλС��
         fontface="italic")              #����������Ϊб��

######��״�ټ����ͼ�״�ټ��صĽ������ö��ڸ�ָ���Ӱ��
result7 <- matrix(0,nrow=20,ncol=6,
              dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                            c("tx1","p1","tx2","p2","tx12","p12")))
ans = NA
#summary(lm(bacs�ܷ�~���������״��ԭ����*��״�ټ���+�Ա�+bmi+����,ms))
for (j in colnames(ms)[43:47]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"*��״�ټ���+�Ա�+bmi+����")
    fit = summary(lm(formula,ms))
    result7[i,c(1,3,5)] = t(fit$coefficients[c(2,3,7),3])
    result7[i,c(2,4,6)] = t(fit$coefficients[c(2,3,7),4])
  }
  ans = cbind(ans,result7)
}
write.csv(ans,"msresult7.csv")

##���ܺ��
ans = NA
summary(lm(bacs�ܷ�~���������״��ԭ����*��״�ټ���+�Ա�+bmi+����,ms))
for (j in colnames(ms)[43:47]){
  results71 <- matrix(0,nrow=5,ncol=3,
                     dimnames=list(colnames(ms)[c(15,30,31,32,37)],
                                   c("tx12","p12")))
  for (i in colnames(ms)[c(15,30,31,32,37)]){
    formula = paste0(i,"~",j,"*��״�ټ���+�Ա�+bmi+����")
    fit = summary(lm(formula,ms))
    results71[i,c(1,2)] = t(fit$coefficients[7,c(3,4)])
  }
  results71 = as.data.frame(results71)
  results71$adjust = p.adjust(results71$p12,method = "bonferroni")
  results71 = as.matrix(results71)
  ans = cbind(ans,results71)
}
write.csv(ans,"msresult7_1.csv")

## WHODAS~TSH/t3/t4*��״�ټ���
result7_2 <- matrix(0,nrow=6,ncol=3,
                    dimnames=list(colnames(ms)[c(9:14)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(9:14)]){
  formula = paste0(i,"~�ټ�״����*��״�ټ���+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
  result7_2[i,"coef"] = fit$coefficients[7,1]
  result7_2[i,"tstat"] = fit$coefficients[7,3]
  result7_2[i,"pvalue"] = fit$coefficients[7,4]
}
result7_2 = as.data.frame(result7_2)
result7_2$adjust = p.adjust(result7_2$pvalue,method = "bonferroni")

write.csv(result7_2,"msresult7_2.csv")

## WHODAS~ft4*��״�ټ���
result7_3 <- matrix(0,nrow=6,ncol=3,
                    dimnames=list(colnames(ms)[c(9:14)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(9:14)]){
  formula = paste0(i,"~�����״����*��״�ټ���+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
  result7_3[i,"coef"] = fit$coefficients[7,1]
  result7_3[i,"tstat"] = fit$coefficients[7,3]
  result7_3[i,"pvalue"] = fit$coefficients[7,4]
}
result7_3 = as.data.frame(result7_3)
result7_3$adjust = p.adjust(result7_3$pvalue,method = "bonferroni")

write.csv(result7_3,"msresult7_3.csv")

## WHODAS~ft4*��״�ټ���
result7_3 <- matrix(0,nrow=6,ncol=3,
                    dimnames=list(colnames(ms)[c(9:14)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(9:14)]){
  formula = paste0(i,"~�����״����*��״�ټ���+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
  result7_3[i,"coef"] = fit$coefficients[7,1]
  result7_3[i,"tstat"] = fit$coefficients[7,3]
  result7_3[i,"pvalue"] = fit$coefficients[7,4]
}
result7_3 = as.data.frame(result7_3)
result7_3$adjust = p.adjust(result7_3$pvalue,method = "bonferroni")

write.csv(result7_3,"msresult7_3.csv")

## WHODAS~T4*��״�ټ���
result7_4 <- matrix(0,nrow=6,ncol=3,
                    dimnames=list(colnames(ms)[c(9:14)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(9:14)]){
  formula = paste0(i,"~��״����*��״�ټ���+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
  result7_4[i,"coef"] = fit$coefficients[7,1]
  result7_4[i,"tstat"] = fit$coefficients[7,3]
  result7_4[i,"pvalue"] = fit$coefficients[7,4]
}
result7_4 = as.data.frame(result7_4)
result7_4$adjust = p.adjust(result7_4$pvalue,method = "bonferroni")

write.csv(result7_4,"msresult7_4.csv")

## PANSS~T4*��״�ټ���
result7_5 <- matrix(0,nrow=3,ncol=3,
                    dimnames=list(colnames(ms)[c(34:36)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(34:36)]){
  formula = paste0(i,"~��״����*��״�ټ���+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
  result7_5[i,"coef"] = fit$coefficients[7,1]
  result7_5[i,"tstat"] = fit$coefficients[7,3]
  result7_5[i,"pvalue"] = fit$coefficients[7,4]
}
result7_5 = as.data.frame(result7_5)
result7_5$adjust = p.adjust(result7_5$pvalue,method = "bonferroni")

write.csv(result7_5,"msresult7_5.csv")

## PANSS~fT4*��״�ټ���
result7_6 <- matrix(0,nrow=3,ncol=3,
                    dimnames=list(colnames(ms)[c(34:36)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(34:36)]){
  formula = paste0(i,"~�����״����*��״�ټ���+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
  result7_6[i,"coef"] = fit$coefficients[7,1]
  result7_6[i,"tstat"] = fit$coefficients[7,3]
  result7_6[i,"pvalue"] = fit$coefficients[7,4]
}
result7_6 = as.data.frame(result7_6)
result7_6$adjust = p.adjust(result7_6$pvalue,method = "bonferroni")

write.csv(result7_6,"msresult7_6.csv")

######����У��BACS����~���������״��ԭ����
result8 <- matrix(0,nrow=7,ncol=1,
              dimnames=list(colnames(ms)[24:30],c("pvalue")))
#summary(lm(bacs�ܷ�~���������״��ԭ����+�Ա�+bmi+����,ms))
for (i in colnames(ms)[24:30]){
  formula = paste0(i,"~���������״��ԭ����+�Ա�+bmi+����")
  fit = summary(lm(formula,ms))
  result8[i,"pvalue"] = fit$coefficients[2,4]
}
result8 = as.data.frame(result8)
result8$adjust = p.adjust(result8$pvalue,method = "bonferroni")
print(result8,digit=3)

######����У��BACS����~PANSS
ans = NA
for (j in colnames(ms)[34:37]){
  result9 <- matrix(0,nrow=7,ncol=1,
                    dimnames=list(colnames(ms)[24:30],c("pvalue")))
  for (i in colnames(ms)[24:30]){
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms))
    result9[i,"pvalue"] = fit$coefficients[2,4]
  }
  result9 = as.data.frame(result9)
  result9$adjust = p.adjust(result9$pvalue,method = "bonferroni")
  ans = cbind(ans,result9)
}
options(digits = 3)
write.csv(ans,"msresult9.csv")




######��ͨ��У��Ļع�����û�����
nPerms = 5000
set.seed(1111)
Th_perm_tr = rlply(nPerms,ms[,24:30])
Th_perm_tr[1]
co_perm1 = matrix(0,nrow=nPerms,ncol=1)
co_perm2 = matrix(0,nrow=nPerms,ncol=1)

summary(lm(bacs�ܷ�~���������״��ԭ����+�Ա�+bmi+����,ms))
summary(lm(bacs����������~���������״��ԭ����+�Ա�+bmi+����,ms))

for (i in 1:nPerms){
  ms1 = Th_perm_tr[[i]]
  fit1 = summary(lm(ms1$bacs�ܷ�~ms$���������״��ԭ����+ms$�Ա�+ms$bmi+ms$����))
  co_perm1[i] = fit1$coefficients[2,1]
  fit2 = summary(lm(ms1$bacs����������~ms$���������״��ԭ����+ms$�Ա�+ms$bmi+ms$����))
  co_perm2[i] = fit2$coefficients[2,1]
}
head(co_perm1)
p_perm11 = sum(abs(co_perm1)>=0.47782)/nPerms
p_perm21 = sum(abs(co_perm2)>=0.54037)/nPerms
print(p_perm11)
print(p_perm21)

######�����Է���

######��ͼBACS����~PANSS
results10 <- matrix(0,nrow=7,ncol=4,
                   dimnames=list(colnames(ms)[24:30],colnames(ms)[34:37]))
#summary(lm(bacs�ܷ�~��ϸ������+bmi+�Ա�+����,ms))
for (i in colnames(ms)[c(24:30)]){
  for (j in colnames(ms)[c(34:37)]){
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms))
    results10[i,j] = fit$coefficients[2,1]
  }
}
#write.csv(results10,"msresult1.csv")
print(results10)
library(ggcorrplot)
ggcorrplot(results10)
#cor��ͼ
painting1 = cor(ms[24:30],ms[34:37], method = "pearson")
print(painting1)
ggcorrplot(painting1)

library(RColorBrewer)
my_color <- brewer.pal(5, "Spectral")# ��ȡ 5 ����ɫ
ggplot(painting1) +
  geom_tile(colour = "black") +
  scale_fill_gradientn(colours = my_color)

######��ͼBACS����~t3
painting2 = cor(ms[24:30],ms$���������״��ԭ����,method = "pearson")
print(painting2)
ggcorrplot(painting2)

######BACS��PANSS�Ĺ�ϵ

###����MMSE����
ms$MMSE = ms$MMSE<24
ms$MMSE = as.factor(ms$MMSE)

###��PANSS,סԺ��Ṧ��,BACS��MMSE�Ĺ�ϵ
results11 <- matrix(0,nrow=12,ncol=3,
                   dimnames=list(colnames(ms)[c(24:30,31,34:37)],
                                 c("coef","tstat","pvalue")))
summary(lm(bacs�ܷ�~MMSE+bmi+�Ա�+����,ms))
for (i in colnames(ms)[c(24:30,31,34:37)]){
  formula = paste0(i,"~MMSE+bmi+�Ա�+����")
  fit = summary(lm(formula,ms))
  results11[i,"coef"] = fit$coefficients["MMSETRUE","Estimate"]
  results11[i,"tstat"] = fit$coefficients["MMSETRUE","t value"]
  results11[i,"pvalue"]  = fit$coefficients["MMSETRUE","Pr(>|t|)"]
}
results11 = as.data.frame(results11)
results11$adjust = p.adjust(results11$pvalue,method = "bonferroni")
#print(results11,digits=4)
write.csv(results11,"msresult11.csv")

###�����Ƿ�մ����飬��BACS��PANSS��סԺ��Ṧ�ܵ�Ӱ��
results12 <- matrix(0,nrow=5,ncol=3,
                  dimnames=list(colnames(ms)[c(31,34:37)],
                                c("coef","tstat","pvalue")))
ans = NA
for (j in colnames(ms)[c(24:30)]){
  results12 <- matrix(0,nrow=5,ncol=3,
                      dimnames=list(colnames(ms)[c(31,34:37)],
                                    c("coef","tstat","pvalue")))
  for (i in colnames(ms)[c(31,34:37)]){
    formula = paste0(i,"~",j,"+�Ա�+bmi+����")
    fit = summary(lm(formula,ms[which(ms$MMSE==FALSE),]))
    results12[i,"coef"] = fit$coefficients[2,1]
    results12[i,"tstat"] = fit$coefficients[2,3]
    results12[i,"pvalue"]  = fit$coefficients[2,4]
  }
  results12 = as.data.frame(results12)
  results12$adjust = p.adjust(results12$pvalue,method = "bonferroni")
  results12 = as.matrix(results12)
  ans = rbind(ans,results12)
}
write.csv(ans,"msresult12.csv")


