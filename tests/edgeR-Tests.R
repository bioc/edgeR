library(edgeR)
options(warnPartialMatchArgs=TRUE,warnPartialMatchAttr=TRUE,warnPartialMatchDollar=TRUE)

set.seed(0); u <- runif(100)

# generate raw counts from NB, create list object
y <- matrix(rnbinom(80,size=5,mu=10),nrow=20)
y <- rbind(0,c(0,0,2,2),y)
rownames(y) <- paste("Tag",1:nrow(y),sep=".")
d <- DGEList(counts=y,group=rep(1:2,each=2),lib.size=1001:1004)

filterByExpr(d)

# estimate common dispersion and find differences in expression
d <- estimateCommonDisp(d)
d$common.dispersion
de <- exactTest(d)
summary(de$table)
topTags(de)

d2 <- estimateTagwiseDisp(d,trend="none",prior.df=20)
summary(d2$tagwise.dispersion)
de <- exactTest(d2,dispersion="common")
topTags(de)

de <- exactTest(d2)
topTags(de)

d2 <- estimateTagwiseDisp(d,trend="movingave",span=0.4,prior.df=20)
summary(d2$tagwise.dispersion)
de <- exactTest(d2)
topTags(de)

summary(exactTest(d2,rejection.region="smallp")$table$PValue)
summary(exactTest(d2,rejection.region="deviance")$table$PValue)

d2 <- estimateTagwiseDisp(d,trend="loess",span=0.8,prior.df=20)
summary(d2$tagwise.dispersion)
de <- exactTest(d2)
topTags(de)

d2 <- estimateTagwiseDisp(d,trend="tricube",span=0.8,prior.df=20)
summary(d2$tagwise.dispersion)
de <- exactTest(d2)
topTags(de)

# mglmOneWay
design <- model.matrix(~group,data=d$samples)
mglmOneWay(d[1:10,],design,dispersion=0.2)
mglmOneWay(d[1:10,],design,dispersion=0)

fit <- glmFit(d,design,dispersion=d$common.dispersion,prior.count=0.5/4)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

fit <- glmFit(d,design,dispersion=d$common.dispersion,prior.count=0.5)
summary(fit$coefficients)

fit <- glmFit(d,design,prior.count=0.5/4)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

dglm <- estimateGLMCommonDisp(d,design)
dglm$common.dispersion
dglm <- estimateGLMTagwiseDisp(dglm,design,prior.df=20)
summary(dglm$tagwise.dispersion)
fit <- glmFit(dglm,design,prior.count=0.5/4)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
dglm <- estimateGLMTrendedDisp(dglm,design)
summary(dglm$trended.dispersion)
dglm <- estimateGLMTrendedDisp(dglm,design,method="power")
summary(dglm$trended.dispersion)
dglm <- estimateGLMTrendedDisp(dglm,design,method="spline")
summary(dglm$trended.dispersion)
dglm <- estimateGLMTrendedDisp(dglm,design,method="bin.spline")
summary(dglm$trended.dispersion)
dglm <- estimateGLMTagwiseDisp(dglm,design,prior.df=20)
summary(dglm$tagwise.dispersion)

dglm2 <- estimateDisp(dglm, design)
summary(dglm2$tagwise.dispersion)
dglm2 <- estimateDisp(dglm, design, prior.df=20)
summary(dglm2$tagwise.dispersion)
dglm2 <- estimateDisp(dglm, design, robust=TRUE)
summary(dglm2$tagwise.dispersion)

# Continuous trend
nlibs <- 3
ntags <- 1000
dispersion.true <- 0.1
# Make first transcript respond to covariate x
x <- 0:2
design <- model.matrix(~x)
beta.true <- cbind(Beta1=2,Beta2=c(2,rep(0,ntags-1)))
mu.true <- 2^(beta.true %*% t(design))
# Generate count data
y <- rnbinom(ntags*nlibs,mu=mu.true,size=1/dispersion.true)
y <- matrix(y,ntags,nlibs)
colnames(y) <- c("x0","x1","x2")
rownames(y) <- paste("Gene",1:ntags,sep="")
d <- DGEList(y)
d <- calcNormFactors(d,method="TMM")
fit <- glmFit(d, design, dispersion=dispersion.true, prior.count=0.5/3)
results <- glmLRT(fit, coef=2)
topTags(results)
d1 <- estimateGLMCommonDisp(d, design, verbose=TRUE)
glmFit(d,design,dispersion=dispersion.true, prior.count=0.5/3)

calcNormFactors(d$counts,method="TMMwsp")

d2 <- estimateDisp(d, design)
summary(d2$tagwise.dispersion)
d2 <- estimateDisp(d, design, prior.df=20)
summary(d2$tagwise.dispersion)
d2 <- estimateDisp(d, design, robust=TRUE)
summary(d2$tagwise.dispersion)

# Exact tests
y <- matrix(rnbinom(20,mu=10,size=3/2),5,4)
group <- factor(c(1,1,2,2))
ys <- splitIntoGroupsPseudo(y,group,pair=c(1,2))
exactTestDoubleTail(ys$y1,ys$y2,dispersion=2/3)

y <- matrix(rnbinom(5*7,mu=10,size=3/2),5,7)
group <- factor(c(1,1,2,2,3,3,3))
ys <- splitIntoGroupsPseudo(y,group,pair=c(1,3))
exactTestDoubleTail(ys$y1,ys$y2,dispersion=2/3)
exactTestBetaApprox(ys$y1,ys$y2,dispersion=2/3)

y[1,3:4] <- 0
design <- model.matrix(~group)
fit <- glmFit(y,design,dispersion=2/3,prior.count=0.5/7)
summary(fit$coefficients)

lrt <- glmLRT(fit,contrast=cbind(c(0,1,0),c(0,0,1)))
topTags(lrt)
design <- model.matrix(~0+group)
fit <- glmFit(y,design,dispersion=2/3,prior.count=0.5/7)
lrt <- glmLRT(fit,contrast=cbind(c(-1,1,0),c(0,-1,1),c(-1,0,1)))
topTags(lrt)

# simple Good-Turing algorithm runs.
test1 <- 1:9
freq1 <- c(2018046, 449721, 188933, 105668, 68379, 48190, 35709, 37710, 22280)
goodTuring(rep(test1, freq1))
test2 <- c(312, 14491, 16401, 65124, 129797, 323321, 366051, 368599, 405261, 604962)
goodTuring(test2)

# Dispersion estimation with fitted values equal to zero
ngenes <- 100
nsamples <- 3
y <- matrix(rnbinom(ngenes*nsamples,size=5,mu=10),ngenes,nsamples)
Group <- factor(c(1,2,2))
design <- model.matrix(~Group)
y[1:5,2:3] <- 0
y <- DGEList(counts=y,group=Group)
y <- estimateCommonDisp(y)
y$common.dispersion
y <- estimateGLMCommonDisp(y,design)
y$common.dispersion
y <- estimateGLMTrendedDisp(y,design)
y$trended.dispersion[1:10]
summary(y$trended.dispersion)
y <- estimateGLMTagwiseDisp(y,design)
y$tagwise.dispersion[1:10]
summary(y$tagwise.dispersion)
y <- estimateDisp(y,design)
y$prior.df
y$common.dispersion
y$trended.dispersion[1:10]
summary(y$trended.dispersion)
y$tagwise.dispersion[1:10]
summary(y$tagwise.dispersion)

# glmQLFit
fit <- glmQLFit(y,design) # legacy=TRUE
fit$dispersion
summary(fit$var.post)
fit$unit.deviance.adj[1:10,]
fit$unit.df.adj[1:10,]
summary(fit$deviance.adj)
summary(fit$df.residual.adj)
