---
layout: post
title: "Developing Distance Based Neighborhood Matrices"
category: R
excerpt_separator: <!--more-->
---

![plot of chunk spantree-plot]({{ site.url }}/images/2016-05-16-distance-based-neighborhoods-spantree-plot-1.png)

For my first few blog posts, I will be porting over some of my previous blogs from another website, and in doing so, updating some of the code base, infrastructure, and content. I figured that I should first start with my blog posts focusing on dealing with autocorrelation utilizing Moran’s Eigenvector Maps (MEM; Dray et al. 2006, Griffith and Peres-Neto 2006).  The objectives of these posts are to introduce the MEM concepts, develop some R code and functions, and to utilize the functions to analyze some data.
<!--more-->

A good place to start is a description of what MEM is and how it is developed.  In short, MEM is multivariate technique that takes a weighted neighborhood matrix and represents it in multivariate space. The advantage of the approach is that it can be incorporated into any linear model based analysis with relative ease. To date, several R packages have been developed recently to perform MEM; however, most lack flexibility and are developed for a certain data type or class. Thus, in these next few posts, I will be developing the functions needed to perform MEM on ordinary matrices and vectors (i.e., data structures that can deal with spatial, temporal, or phylogenetic information).  In addition, I will be developing some methods to select eigenvectors to include in linear models.

This first post will focus on taking a distance matrix and creating a neighborhood matrix from it. This neighborhood matrix can than be used to:

1. Test for autocorrelation using Moran’s I
2. Develop spatial, temporal, spatiotemporal, or phylogenetic eigenvectors that represent autocorrelation patterns.

There are a variety of rules that we can use to define these neighborhood matrices, especially if we are dealing with lattice data (i.e., data located in regions). For these cases, we can use rules based on chess piece movements (e.g., rook, queen) to define regions (e.g., counties) that are neighboring. However, often we deal with data within a domain where our data does not explicitly fit into regions, but rather they are points in space, time or even from an evolutionary spectrum. For these cases, we need some way of defining whether a site, time, or species is a neighbor or not with another site, time, or species and a way of expressing the difficultly of exchange between other sites, times, or species. This is where the weighted neighborhood matrix comes in, and this will be the topic of this post.

If we have data that are points in space or time, we can extract a distance matrix from them (this also can be done for phylogenetic distance with the proper genetic information). For example if we have two vectors representing an $X$ and $Y$ coordinate (on a flat plane [if $X$ and $Y$ are latitude and longitude we might need to use a function like geoXY()]), we can calculate a distance matrix in R using the `dist` function:


```r
dist(cbind(x,y))
```

Here we have combined the two vectors into a matrix using the `cbind` function. By default we receive the lower triangle of the distance matrix, but we can see the entire matrix using the command:


```r
dist(cbind(x,y),upper=T,diag=T)
```

or by changing the class of the object created


```r
as.matrix(dist(cbind(x,y)))
```

If this distance matrix ($\mathbf{D}$) is stored as an object in `R`, we can use a variety of rules to calculate a weighted neighborhood matrix. Below, I focus on four different functions that describe the increase in difficulty of exchange between points or species. These functions are the ones described by Dray et al. (2006) and include: distance-based Moran's Eigenvector Maps (dbmem), Linear Function, Concave-Down, and Concave-Up.

For all these methods, we need to derive two separate matrices: (1) a connectivity matrix and (2) an edge weight matrix. The connectivity matrix ($\mathbf{B}$) simply contains information on whether points are connected and is binary, with 1s representing connection and 0s representing unconnected sites. The edge weighting matrix ($\mathbf{A}$) is a little more complex and represents the difficulty of exchange between all points in the sample domain.

For the connectivity matrix, all we need is a threshold in which to say sites are no longer connected. Of course we want our function to be flexible, so we should allow this to be a user defined threshold, but there is a standard threshold used for dbmem. This is usually defined as the maximum distance on a minimum span tree. This span tree is simply a tree that connects all sites together with edges, and the minimum span tree is the tree with the lowest cost (i.e., represents the least difficult exchange). Luckily, the `R` package `vegan` provides an easy way of calculating the minimum span tree, of which we only need to find the maximum.


```r
library(vegan)
thresh=max(spantree(D)$dist)
```

Once the threshold is defined, we can simply calculate \mathbf{B} using an `apply` statement:


```r
B<-apply(D,2,function(x)ifelse(x>thresh,0,1))
```

Calculating $\mathbf{A}$ can be done similarly, but the weighting depends on which method is used. For dbmem we need to define $\mathbf{A}$ such that:

$$\mathbf{A}=[a_{ij}]=1-\Big (\frac{d_{ij}}{4 \times thresh} \Big )^2, $$

where $d_{ij}$ indexes the elements of $\mathbf{D}$. The linear function is defined similarly, but no threshold is specified:

$$\mathbf{A}=[a_{ij}]=1-\Big (\frac{d_{ij}}{max(\mathbf{D})} \Big ).$$

Likewise, the concave-down function is similar to the linear function, with the latter half raised to the power of alpha:

$$\mathbf{A}=[a_{ij}]=1-\Big (\frac{d_{ij}}{max(\mathbf{D})} \Big )^\alpha.$$ 

Finally, the concave-up function is defined as:

$$\mathbf{A}=[a_{ij}]=1/d_{ij}^\beta.$$

The concave-up function represents a similarity that decreases less rapidly for larger values of $\mathbf{D}$. The methods can be easily calculated with an apply function in `R` (see code for full function below).

Finally, to compute the weighted neighborhood we just need to take the Hadamard product of $mathbf{B}$ and $\mathbf{A}$. If `R` this is simply,


```r
W<-B*A
```

However, I should note that the diagonal of $\mathbf{W}$ will be ones representing a connection. This is good if we are going to try to derive more complex relationships with other weight matrices (e.g., spatiotemporal designs), but for the most part we would like to have zeros on the diagonal. This can be done simply with the statement:


```r
diag(W)<-0
```

Below is a script that is available on github where I have combined all of these weight functions into a general `R` function named `weightmatrix`. This function requires a distance matrix, the method to be specified, and a few other arguments for flexibility of the function. By default the function performs dbmem weighting and uses a threshold based on a minimum span tree.



```r
weightmatrix<-function (D, method = c("dbmem", "linear", "concave-up", "concave-down","connectivity"), thresh = NULL, alpha = NULL, beta = NULL, diag.mat = 0) {
  if (class(D) == "dist")
    D = as.matrix(D)
  require(vegan)
  n <- dim(D)[1]
  if (is.null(thresh))
    thresh = max(spantree(D)$dist)
  out = A = (B <- array(dim = dim(D)))
  if (length(method) > 1)
    method = "dbmem"
  B<-apply(D,2,function(x)ifelse(x>thresh,0,1))
  if(method=="dbmem") A<-apply(D,2,function(x)1 - ((x/(4 * thresh))^2))
  else if(method=="linear")A<-apply(D,2,function(x)1 - (x/max(D)))
  else if (method == "concave-down") {
    if (is.null(alpha)) stop("alpha must be specified for the method 'concave-down'")
    A<-apply(D,2,function(x)1 - (x/max(D))^alpha)
  }
  else if (method == "concave-up") {
    if (is.null(beta)) stop("beta must be specified for the method 'concave-up")
    A<-apply(D,2,function(x)1/x^beta)
  }
  else if (method == "connectivity")A<-B
  else stop("method must be 'dbmem', 'linear', 'concave-up', 'concave-down', or 'connectivity'")
  out <- A * B
  if (diag.mat == 0) {
    diag(out) <- 0
  }
  return(out)
}
```

Now, let’s try the function. We need to simulate some data and run the function. For sake of simplicity, we will perform the `dbmem` weighting.


```r
D<-as.matrix(dist(cbind(runif(100,0,100),runif(100,0,100))))
W<-weightmatrix(D,method="dbmem")
W
```

Everything seems to work!

As with most of the functions created on this blog, the weightmatrix function is available in my github repository and can be installed (along with other functions from my blog) using the following `R` commands:


```r
install.packages("devtools")
library(devtools)
install_github("quanteco/quanteco")
library(quanteco)
```

The benefit of installing the package, instead of coping the code, is that you will also get the associated help files and examples.

# References

Dray, S., P. Legendre, and P. R. Peres-Neto. 2006. Spatial modeling: a comprehensive framework for principal coordinate analysis of neighborhood matrices (PCNM). Ecological Modeling 196:483-493.

Griffith, D. A., and P. R. Peres-Neto. 2006. Spatial modeling in ecology: the flexibility of eigenfunction spatial analyses. Ecology 87:2603–2613.
