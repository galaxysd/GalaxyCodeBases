#!/usr/bin/env littler

#install.packages('e1071')
require('e1071')
#mdat <- matrix(c(21,.1,.11,.12,.14, 23,.15,.16,.17,.15, 27,.2,.19,.18,.3, 32,.3,.33,.32,.34, 36,.4,.41,.42,.4, 42,.5,.51,.52,.54), 5, 5, TRUE, dimnames = list(c(),c("Age", "Pos1", "Pos2","Pos3","Pos4")))
mdat <- read.delim('dosvm.txt',comment.char='#')
x <- subset(mdat, select = -Age)
y <- mdat[,1]
model <- svm(Age ~ ., data = mdat, type='eps-regression', gamma=0.1, cost=2, epsilon = 0.1)
model2 <- svm(Age ~ a11044875 + a11044880, data = mdat, type='eps-regression', gamma=0.1, cost=2, epsilon = 0.1)
model3 <- svm(Age ~ a11044875 + a11044880 + b11044877 + b11044888 + b11044894, data = mdat, type='eps-regression', gamma=0.1, cost=2, epsilon = 0.1)
pred_result <- predict(model, x)
pred_result2 <- predict(model2, x)
pred_result3 <- predict(model3, x)

#print(table(pred_result,y))

print(model)
print(y)
print(pred_result)
print(pred_result2)
print(pred_result3)

#plot(model, mdat, Pos1 ~ Pos2,color.palette = terrain.colors)
