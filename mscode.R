library(plyr)
library(dplyr) #调用select
library(psych) #相关系数检验
library(pacman) #缺失值处理
library(mice)
library("VIM")
library(pheatmap)
library(car)
library(ggplot2)
library(corrplot)
library(RColorBrewer)

######读入数据

setwd("d:/fx")
ms = read.table(file = "ms.csv", header = T , 
                     sep = "," , fill = TRUE , encoding = "UTF-8")
colnames(ms)[1] = "病程"
head(ms)
ms = ms[1:45,]

######修改BACS相关数据为zscore
ms$bacs语言记忆 = scale(ms$bacs语言记忆)
ms$bacs数字序列 = scale(ms$bacs数字序列)
ms$bacs代币动作 = scale(ms$bacs代币动作)
ms$bacs伦敦塔 = scale(ms$bacs伦敦塔)
ms$bacs符号替代 = scale(ms$bacs符号替代)
ms$bacs语义流畅度 = scale(ms$bacs语义流畅度)
ms$bacs总分 = scale(ms$bacs总分)

######BACS与慢性病、PANSS、住院社会功能、血清指标的相关性
results1 <- matrix(0,nrow=48,ncol=7,
                   dimnames=list(rep(c("t","p"),24),(colnames(ms)[c(24:30)])))
#summary(lm(bacs总分~白细胞计数+bmi+性别+病程,ms))
for (i in colnames(ms)[c(24:30)]){
  s = 0
  for (j in colnames(ms)[c(18:23,34:37,31,41:53)]){
    s = s+1
    formula = paste0(i,"~",j,"+性别+bmi+病程")
    fit = summary(lm(formula,ms))
    results1[(s-1)*2+1,i] = fit$coefficients[2,3]
    results1[(s-1)*2+2,i] = fit$coefficients[2,4]
  }
}
write.csv(results1,"msresult1.csv")


######住院社会功能与慢性病、panss、血清的相关性
results2 <- matrix(0,nrow=16,ncol=2,
                   dimnames=list(colnames(ms)[c(18:23,34:37,48:53)],c("t","p")))
for (i in colnames(ms)[c(18:23,34:37,48:53)]){
  formula = paste0("住院社会功能~",i,"+性别+bmi+病程")
  fit = summary(lm(formula,ms))
    results2[i,1] = fit$coefficients[2,3]
    results2[i,2] = fit$coefficients[2,4]
}
write.csv(results2,"msresult2.csv")

######CRP与BACS、慢性病、PANSS、住院社会功能、血清指标的相关性
results3 <- matrix(0,nrow=25,ncol=2,
                   dimnames=list(colnames(ms)[c(18:32,34:37,48:53)],
                                 c("t","p")))
for (i in colnames(ms)[c(18:32,34:37,48:53)]){
  formula = paste0("crp~",i,"+性别+bmi+病程")
  fit = summary(lm(formula,ms))
  results3[i,1] = fit$coefficients[2,3]
  results3[i,2] = fit$coefficients[2,4]
}
write.csv(results3,"msresult3.csv")

######是否患有甲状腺疾病对于WHOQOL、住院生活质量、WHO-DAS、BACS及PANSS的影响
ms$甲状腺疾病 = as.factor(ms$甲状腺疾病)
result3_5 <- matrix(0,nrow=20,ncol=3,
                  dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],c("coef","tstat","pvalue")))

summary(lm(bacs总分~甲状腺疾病+bmi+性别+病程,ms))
for (i in colnames(ms)[c(9:15,24:32,34:37)]){
  formula = paste0(i,"~甲状腺疾病+性别+bmi+病程")
  fit = summary(lm(formula,ms))
  result3_5[i,"coef"] = fit$coefficients[2,1]
  result3_5[i,"tstat"] = fit$coefficients[2,3]
  result3_5[i,"pvalue"] = fit$coefficients[2,4]
}

write.csv(result3_5,"msresult3_5.csv")



######各甲状腺激素对于WHOQOL、住院生活质量、WHO-DAS、BACS及PANSS的影响
result4 <- matrix(0,nrow=20,ncol=3,
                     dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                                   c("coef","tstat","pvalue")))
ans = NA
#summary(lm(bacs总分~游离甲状腺素+bmi+性别+病程,ms))
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+性别+bmi+病程")
    fit = summary(lm(formula,ms))
    result4[i,"coef"] = fit$coefficients[2,1]
    result4[i,"tstat"] = fit$coefficients[2,3]
    result4[i,"pvalue"]  = fit$coefficients[2,4]
  }
  ans = cbind(ans,result4)
}
write.csv(ans,"msresult4.csv")

#result4相关系数绘图
vis4 = matrix(0,nrow=20,ncol=5,
              dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                            colnames(ms)[c(43:47)]))
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+性别+bmi+病程")
    fit = summary(lm(formula,ms))
    vis4[i,j] = fit$coefficients[2,4]
  }
}
# 将单元格的高度和宽度都设置为20
pheatmap(-log(vis4,base=10), border_color = "black", #边框线为黑色
         display_numbers = F,         #热图格子中显示相应的数值
         number_color = "black",         #字体颜色为黑色
         fontsize=10,                    #字体大小为10
         cluster_rows=F,                 #不进行聚类
         cluster_cols = F,
         number_format = "%.3f",         #保留一位小数
         fontface="italic")              #将字体设置为斜体

######(患有甲状腺疾病)各甲状腺激素对于WHOQOL、住院生活质量、WHO-DAS、BACS及PANSS的影响
result5 <- matrix(0,nrow=20,ncol=3,
                  dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                                c("coef","tstat","pvalue")))
ans = NA
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+性别+bmi+病程")
    fit = summary(lm(formula,ms[which(ms$甲状腺疾病==1),]))
    result5[i,"coef"] = fit$coefficients[2,1]
    result5[i,"tstat"] = fit$coefficients[2,3]
    result5[i,"pvalue"]  = fit$coefficients[2,4]
  }
  ans = cbind(ans,result5)
}
write.csv(ans,"msresult5.csv")

#result5相关系数绘图
vis5 = matrix(0,nrow=20,ncol=5,
              dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                            colnames(ms)[c(43:47)]))
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+性别+bmi+病程")
    fit = summary(lm(formula,ms[which(ms$甲状腺疾病==1),]))
    vis5[i,j] = fit$coefficients[2,4]
  }
}
# 将单元格的高度和宽度都设置为20
pheatmap(-log(vis5,base=10), border_color = "black", #边框线为黑色
         display_numbers = F,         #热图格子中显示相应的数值
         number_color = "black",         #字体颜色为黑色
         fontsize=10,                    #字体大小为10
         cluster_rows=F,                 #不进行聚类
         cluster_cols = F,
         number_format = "%.3f",         #保留一位小数
         fontface="italic")              #将字体设置为斜体

######(未患有甲状腺疾病)各甲状腺激素对于WHOQOL、住院生活质量、WHO-DAS、BACS及PANSS的影响
result6 <- matrix(0,nrow=20,ncol=3,
                  dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                                c("coef","tstat","pvalue")))
ans = NA
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+性别+bmi+病程")
    fit = summary(lm(formula,ms[which(ms$甲状腺疾病==0),]))
    result6[i,"coef"] = fit$coefficients[2,1]
    result6[i,"tstat"] = fit$coefficients[2,3]
    result6[i,"pvalue"]  = fit$coefficients[2,4]
  }
  ans = cbind(ans,result6)
}
write.csv(ans,"msresult6.csv")
#result6相关系数绘图
vis6 = matrix(0,nrow=20,ncol=5,
              dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                            colnames(ms)[c(43:47)]))
for (j in colnames(ms)[c(43:47)]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"+性别+bmi+病程")
    fit = summary(lm(formula,ms[which(ms$甲状腺疾病==0),]))
    vis6[i,j] = fit$coefficients[2,4]
  }
}
# 将单元格的高度和宽度都设置为20
pheatmap(-log(vis6,base=10), border_color = "black", #边框线为黑色
         display_numbers = F,         #热图格子中显示相应的数值
         number_color = "black",         #字体颜色为黑色
         fontsize=10,                    #字体大小为10
         cluster_rows=F,                 #不进行聚类
         cluster_cols = F,
         number_format = "%.3f",         #保留一位小数
         fontface="italic")              #将字体设置为斜体

######甲状腺疾病和甲状腺激素的交互作用对于各指标的影响
result7 <- matrix(0,nrow=20,ncol=6,
              dimnames=list(colnames(ms)[c(9:15,24:32,34:37)],
                            c("tx1","p1","tx2","p2","tx12","p12")))
ans = NA
#summary(lm(bacs总分~游离三碘甲状腺原氨酸*甲状腺疾病+性别+bmi+病程,ms))
for (j in colnames(ms)[43:47]){
  for (i in colnames(ms)[c(9:15,24:32,34:37)]){
    formula = paste0(i,"~",j,"*甲状腺疾病+性别+bmi+病程")
    fit = summary(lm(formula,ms))
    result7[i,c(1,3,5)] = t(fit$coefficients[c(2,3,7),3])
    result7[i,c(2,4,6)] = t(fit$coefficients[c(2,3,7),4])
  }
  ans = cbind(ans,result7)
}
write.csv(ans,"msresult7.csv")

##先总后分
ans = NA
summary(lm(bacs总分~游离三碘甲状腺原氨酸*甲状腺疾病+性别+bmi+病程,ms))
for (j in colnames(ms)[43:47]){
  results71 <- matrix(0,nrow=5,ncol=3,
                     dimnames=list(colnames(ms)[c(15,30,31,32,37)],
                                   c("tx12","p12")))
  for (i in colnames(ms)[c(15,30,31,32,37)]){
    formula = paste0(i,"~",j,"*甲状腺疾病+性别+bmi+病程")
    fit = summary(lm(formula,ms))
    results71[i,c(1,2)] = t(fit$coefficients[7,c(3,4)])
  }
  results71 = as.data.frame(results71)
  results71$adjust = p.adjust(results71$p12,method = "bonferroni")
  results71 = as.matrix(results71)
  ans = cbind(ans,results71)
}
write.csv(ans,"msresult7_1.csv")

## WHODAS~TSH/t3/t4*甲状腺疾病
result7_2 <- matrix(0,nrow=6,ncol=3,
                    dimnames=list(colnames(ms)[c(9:14)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(9:14)]){
  formula = paste0(i,"~促甲状腺素*甲状腺疾病+性别+bmi+病程")
  fit = summary(lm(formula,ms))
  result7_2[i,"coef"] = fit$coefficients[7,1]
  result7_2[i,"tstat"] = fit$coefficients[7,3]
  result7_2[i,"pvalue"] = fit$coefficients[7,4]
}
result7_2 = as.data.frame(result7_2)
result7_2$adjust = p.adjust(result7_2$pvalue,method = "bonferroni")

write.csv(result7_2,"msresult7_2.csv")

## WHODAS~ft4*甲状腺疾病
result7_3 <- matrix(0,nrow=6,ncol=3,
                    dimnames=list(colnames(ms)[c(9:14)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(9:14)]){
  formula = paste0(i,"~游离甲状腺素*甲状腺疾病+性别+bmi+病程")
  fit = summary(lm(formula,ms))
  result7_3[i,"coef"] = fit$coefficients[7,1]
  result7_3[i,"tstat"] = fit$coefficients[7,3]
  result7_3[i,"pvalue"] = fit$coefficients[7,4]
}
result7_3 = as.data.frame(result7_3)
result7_3$adjust = p.adjust(result7_3$pvalue,method = "bonferroni")

write.csv(result7_3,"msresult7_3.csv")

## WHODAS~ft4*甲状腺疾病
result7_3 <- matrix(0,nrow=6,ncol=3,
                    dimnames=list(colnames(ms)[c(9:14)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(9:14)]){
  formula = paste0(i,"~游离甲状腺素*甲状腺疾病+性别+bmi+病程")
  fit = summary(lm(formula,ms))
  result7_3[i,"coef"] = fit$coefficients[7,1]
  result7_3[i,"tstat"] = fit$coefficients[7,3]
  result7_3[i,"pvalue"] = fit$coefficients[7,4]
}
result7_3 = as.data.frame(result7_3)
result7_3$adjust = p.adjust(result7_3$pvalue,method = "bonferroni")

write.csv(result7_3,"msresult7_3.csv")

## WHODAS~T4*甲状腺疾病
result7_4 <- matrix(0,nrow=6,ncol=3,
                    dimnames=list(colnames(ms)[c(9:14)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(9:14)]){
  formula = paste0(i,"~甲状腺素*甲状腺疾病+性别+bmi+病程")
  fit = summary(lm(formula,ms))
  result7_4[i,"coef"] = fit$coefficients[7,1]
  result7_4[i,"tstat"] = fit$coefficients[7,3]
  result7_4[i,"pvalue"] = fit$coefficients[7,4]
}
result7_4 = as.data.frame(result7_4)
result7_4$adjust = p.adjust(result7_4$pvalue,method = "bonferroni")

write.csv(result7_4,"msresult7_4.csv")

## PANSS~T4*甲状腺疾病
result7_5 <- matrix(0,nrow=3,ncol=3,
                    dimnames=list(colnames(ms)[c(34:36)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(34:36)]){
  formula = paste0(i,"~甲状腺素*甲状腺疾病+性别+bmi+病程")
  fit = summary(lm(formula,ms))
  result7_5[i,"coef"] = fit$coefficients[7,1]
  result7_5[i,"tstat"] = fit$coefficients[7,3]
  result7_5[i,"pvalue"] = fit$coefficients[7,4]
}
result7_5 = as.data.frame(result7_5)
result7_5$adjust = p.adjust(result7_5$pvalue,method = "bonferroni")

write.csv(result7_5,"msresult7_5.csv")

## PANSS~fT4*甲状腺疾病
result7_6 <- matrix(0,nrow=3,ncol=3,
                    dimnames=list(colnames(ms)[c(34:36)],c("coef","tstat","pvalue")))

for (i in colnames(ms)[c(34:36)]){
  formula = paste0(i,"~游离甲状腺素*甲状腺疾病+性别+bmi+病程")
  fit = summary(lm(formula,ms))
  result7_6[i,"coef"] = fit$coefficients[7,1]
  result7_6[i,"tstat"] = fit$coefficients[7,3]
  result7_6[i,"pvalue"] = fit$coefficients[7,4]
}
result7_6 = as.data.frame(result7_6)
result7_6$adjust = p.adjust(result7_6$pvalue,method = "bonferroni")

write.csv(result7_6,"msresult7_6.csv")

######分组校验BACS各项~游离三碘甲状腺原氨酸
result8 <- matrix(0,nrow=7,ncol=1,
              dimnames=list(colnames(ms)[24:30],c("pvalue")))
#summary(lm(bacs总分~游离三碘甲状腺原氨酸+性别+bmi+病程,ms))
for (i in colnames(ms)[24:30]){
  formula = paste0(i,"~游离三碘甲状腺原氨酸+性别+bmi+病程")
  fit = summary(lm(formula,ms))
  result8[i,"pvalue"] = fit$coefficients[2,4]
}
result8 = as.data.frame(result8)
result8$adjust = p.adjust(result8$pvalue,method = "bonferroni")
print(result8,digit=3)

######分组校验BACS各项~PANSS
ans = NA
for (j in colnames(ms)[34:37]){
  result9 <- matrix(0,nrow=7,ncol=1,
                    dimnames=list(colnames(ms)[24:30],c("pvalue")))
  for (i in colnames(ms)[24:30]){
    formula = paste0(i,"~",j,"+性别+bmi+病程")
    fit = summary(lm(formula,ms))
    result9[i,"pvalue"] = fit$coefficients[2,4]
  }
  result9 = as.data.frame(result9)
  result9$adjust = p.adjust(result9$pvalue,method = "bonferroni")
  ans = cbind(ans,result9)
}
options(digits = 3)
write.csv(ans,"msresult9.csv")




######对通过校验的回归进行置换检验
nPerms = 5000
set.seed(1111)
Th_perm_tr = rlply(nPerms,ms[,24:30])
Th_perm_tr[1]
co_perm1 = matrix(0,nrow=nPerms,ncol=1)
co_perm2 = matrix(0,nrow=nPerms,ncol=1)

summary(lm(bacs总分~游离三碘甲状腺原氨酸+性别+bmi+病程,ms))
summary(lm(bacs语义流畅度~游离三碘甲状腺原氨酸+性别+bmi+病程,ms))

for (i in 1:nPerms){
  ms1 = Th_perm_tr[[i]]
  fit1 = summary(lm(ms1$bacs总分~ms$游离三碘甲状腺原氨酸+ms$性别+ms$bmi+ms$病程))
  co_perm1[i] = fit1$coefficients[2,1]
  fit2 = summary(lm(ms1$bacs语义流畅度~ms$游离三碘甲状腺原氨酸+ms$性别+ms$bmi+ms$病程))
  co_perm2[i] = fit2$coefficients[2,1]
}
head(co_perm1)
p_perm11 = sum(abs(co_perm1)>=0.47782)/nPerms
p_perm21 = sum(abs(co_perm2)>=0.54037)/nPerms
print(p_perm11)
print(p_perm21)

######敏感性分析

######绘图BACS各项~PANSS
results10 <- matrix(0,nrow=7,ncol=4,
                   dimnames=list(colnames(ms)[24:30],colnames(ms)[34:37]))
#summary(lm(bacs总分~白细胞计数+bmi+性别+病程,ms))
for (i in colnames(ms)[c(24:30)]){
  for (j in colnames(ms)[c(34:37)]){
    formula = paste0(i,"~",j,"+性别+bmi+病程")
    fit = summary(lm(formula,ms))
    results10[i,j] = fit$coefficients[2,1]
  }
}
#write.csv(results10,"msresult1.csv")
print(results10)
library(ggcorrplot)
ggcorrplot(results10)
#cor绘图
painting1 = cor(ms[24:30],ms[34:37], method = "pearson")
print(painting1)
ggcorrplot(painting1)

library(RColorBrewer)
my_color <- brewer.pal(5, "Spectral")# 获取 5 个颜色
ggplot(painting1) +
  geom_tile(colour = "black") +
  scale_fill_gradientn(colours = my_color)

######绘图BACS各项~t3
painting2 = cor(ms[24:30],ms$游离三碘甲状腺原氨酸,method = "pearson")
print(painting2)
ggcorrplot(painting2)

######BACS和PANSS的关系

###按照MMSE分组
ms$MMSE = ms$MMSE<24
ms$MMSE = as.factor(ms$MMSE)

###看PANSS,住院社会功能,BACS和MMSE的关系
results11 <- matrix(0,nrow=12,ncol=3,
                   dimnames=list(colnames(ms)[c(24:30,31,34:37)],
                                 c("coef","tstat","pvalue")))
summary(lm(bacs总分~MMSE+bmi+性别+病程,ms))
for (i in colnames(ms)[c(24:30,31,34:37)]){
  formula = paste0(i,"~MMSE+bmi+性别+病程")
  fit = summary(lm(formula,ms))
  results11[i,"coef"] = fit$coefficients["MMSETRUE","Estimate"]
  results11[i,"tstat"] = fit$coefficients["MMSETRUE","t value"]
  results11[i,"pvalue"]  = fit$coefficients["MMSETRUE","Pr(>|t|)"]
}
results11 = as.data.frame(results11)
results11$adjust = p.adjust(results11$pvalue,method = "bonferroni")
#print(results11,digits=4)
write.csv(results11,"msresult11.csv")

###按照是否痴呆分组，看BACS对PANSS和住院社会功能的影响
results12 <- matrix(0,nrow=5,ncol=3,
                  dimnames=list(colnames(ms)[c(31,34:37)],
                                c("coef","tstat","pvalue")))
ans = NA
for (j in colnames(ms)[c(24:30)]){
  results12 <- matrix(0,nrow=5,ncol=3,
                      dimnames=list(colnames(ms)[c(31,34:37)],
                                    c("coef","tstat","pvalue")))
  for (i in colnames(ms)[c(31,34:37)]){
    formula = paste0(i,"~",j,"+性别+bmi+病程")
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



