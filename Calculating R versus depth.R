rm(list=ls())
require(pracma)

rH=function(points,h){
	x=points[,1]
	y=points[,2]
	# To find left points above and below water height
	leftu=min(which(y <= h))-1
	leftb=min(which(y <= h))
	# To find right points above and below water height
	rightu=max(which(y <= h))+1
	rightb=max(which(y <= h))
	#Make a table
	xy.table=cbind(x[leftu:rightu],y[leftu:rightu])
	#Name columns 
	colnames(xy.table)=c("x","y")
	#Get last data point
	lastone=length(xy.table[,1])
	
	#Get x and y for last two data points
	x2=xy.table[lastone,1]
	y2=xy.table[lastone,2]
	
	x1=xy.table[lastone-1,1]
	y1=xy.table[lastone-2,2]
	
	#Equation for right side of cross section (xright,h)
	x.right=x1+(x2-x1)*((h-y1)/(y2-y1))
	y.right=h
	
	xy.table[lastone,1]=x.right
	xy.table[lastone,2]=y.right
	
	#Equation for left side of cross section (xleft,h)
	x1=xy.table[1,1]
	y1=xy.table[1,2]
	
	x2=xy.table[2,1]
	y2=xy.table[2,2]
	
	x.left=x1+(x2-x1)*((h-y1)/(y2-y1))
	#x.left
	y.left=h
	
	xy.table[1,1]=x.left
	xy.table[1,2]=y.left
	
	#plot(xy.table, type="l")
	#print(xy.table)

	perimeter = perimeter(xy.table)	
	area = polyArea(xy.table)
	return(area/perimeter)
}


perimeter = function(points){
	x=points[,1]
	y=points[,2]
	perimetr= 0
	for(i in 1:(length(x)-1)){
	perimetr= perimetr + pyth(x[i],y[i],x[i+1],y[i+1])
	}
	return(perimetr)
}

pyth=function(x1,y1,x2,y2){
	return(sqrt((x1-x2)^2+(y1-y2)^2))
}

polyAreaH <- function(points,h){
	x <- points[,1]
	y <- points[,2]
	area = 0
	temp1 <- which(y <= h)
	left <- min(temp1)
	right <- max(temp1)
	if(x[left-1]-x[left]==0){
		leftx <- x[left] } else {
		leftx <- x[left]+(h-y[left])/(y[left-1]-y[left])*(x[left-1]-x[left])
	}
	if(x[right+1]-x[right]==0){
		rightx <- x[right] } else {
		rightx <- x[right]+(h-y[right])/(y[right+1]-y[right])*(x[right+1]-x[right])
	}
	xh <- c(leftx,x[left:right],rightx)
	yh <- c(h,y[left:right],h)
	xyh <- cbind(xh,yh)
	return(c(polyArea(xyh),polyPeri(xyh)))  
}

polyArea <- function(points){
	x <- points[,1]
	y <- points[,2]
	area = 0
	j <- length(x)
	for(i in 1:j){
		area <- area + (x[j] + x[i]) * (y[j] - y[i])
		j = i
	}
	return(abs(area/2))
}

# --------------- End Functions -----------------

# --------------- Begin Data ------------------
xs_data = read.csv("xs_1-5_data.csv",header=T)
long_data = read.csv("xs_1-5_properties.csv",header=T)
xs_depths = read.csv("probe_depths_april.csv",header=T)

attach(xs_data)

# --------------- End Data ------------------

# Load XS Data ----------------------------------
#Cross Section at Probe 5
xs_data5 = cbind(p5_x, p5_y)
xs_data5 = xs_data5[c(-27:-47),]

#Cross Section at Probe 4
xs_data4 = cbind(p4_x,p4_y)
xs_data4 = xs_data4[c(-25:-47),]

#Cross Section at Probe 3
xs_data3 = cbind(p3_x,p3_y)
xs_data3 = xs_data3[c(-35:-47),]

#Cross Section at Probe 2
xs_data2 = cbind(p2_x,p2_y)
xs_data2 = xs_data2[c(-27:-47),]

#Cross Section at Probe 1
xs_data1 = cbind(p1_x,p1_y)
xs_data1 = xs_data1[c(-39:-47),]

f = list(xs_data5,xs_data4,xs_data3,xs_data2,xs_data1)

# End Load XS Data ------------------------------

# --------------- Parameters ------------------

# --------------- Calculations ------------------

xs.rh <- matrix(NA,nrow=50,ncol=6)
depths <- seq(0.1,9.9,0.2)
length(depths)
xs.rh[,1] = depths

for(i in 1:5){
	xs_xy = as.data.frame(f[i])
		for(j in 1:50){
			 xs.rh[j,i + 1]<- rH(xs_xy,xs.rh[j,1])
			 print(xs.rh[j,1])
		}
}

xs.rh

plot(xs.rh[,1],xs.rh[,6])