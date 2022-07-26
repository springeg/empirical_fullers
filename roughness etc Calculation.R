rm(list=ls())
require(pracma)

a.p=function(points,h){
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
	x.left
	y.left=h
	
	xy.table[1,1]=x.left
	xy.table[1,2]=y.left
	
	#plot(xy.table, type="l")
	#print(xy.table)

	perimeter = perimeter(xy.table)	
	return(c(polyArea(xy.table),perimeter))
	#return(c(polyArea(xy.table[,1],xy.table[,2]),perimeter))
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

solver=function(Q, H1, points, z2, n, l1, l2){
	it=NULL # Want to be zero, is difference from zero
	guess=0.1	
	
	for(i in 1:10000){
		#print(guess)
		ap=a.p(points,guess)
		A2=ap[1]
		u=Q/A2 #velocity
		r=ap[1]/ap[2] #hydraulic radius
		it[i]=H1-((Q/A2)^2/(2*9.81))-guess-z2+(l2-l1)*(n^2*u^2/r^(4/3))
		#print(c(guess,it[i]))

		if(abs(it[i])< 0.01){
			return(guess)
			}
		guess=guess+0.001
		}

return(it)
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
# Discharge
Qq <- c(0.48,0.86,1.29,1.66,0.82,0.45)

# Columns
d.run = 4
# Discharge
Q <- Qq[d.run]
# Probes
probes <- c(5:1)
# Elevations
z <- long_data[,2]
# True depths at XSs
true_depths <- xs_depths[,d.run + 1]
xx <- long_data[,1]
wb <- c(0.98,1.85,0.7,1.06,1.48)
area <- NULL
peri <- NULL
uu <- NULL
h.ead <- NULL
nn <- NULL
so <- NULL
sf = NULL
sws= NULL

# --------------- Calculations ------------------

for(i in 1:4){
	so[i] = (z[i+1]-z[i])/(xx[i+1]-xx[i])
}
print(c("so"))
print(c(round(so,3)))


for(i in 1:5){
	xs_xy = as.data.frame(f[i])
	ap <- a.p(xs_xy,true_depths[i])
	area[i] <- ap[1]
	peri[i] <- ap[2]
	uu[i] <- Q / ap[1]	
}

# Hydraulic radius
rh <- area/peri

# Total head
h.ead <- z + uu^2/(2*9.81) + true_depths
WSE <- z + true_depths

h.ave <- rep(NA,4)
for(i in 1:4){
	h.ave[i] <- mean(c(true_depths[i],true_depths[i+1]))
}
print(c("h.ave"))
print(c(round(h.ave,2)))

wb.ave <- rep(NA,4)
for(i in 1:4){
	wb.ave[i] <- mean(c(wb[i],wb[i+1]))
}
print(c("w.ave"))
print(c(round(w.ave,2)))

h.ave <- rep(NA,4)
for(i in 1:4){
	h.ave[i] <- mean(c(true_depths[i],true_depths[i+1]))
}
print(c("h.ave"))
print(c(round(h.ave,2)))

uu.ave <- rep(NA,4)
for(i in 1:4){
	uu.ave[i] <- mean(c(uu[i],uu[i+1]))
}
print(c("uu.ave"))
print(c(round(uu.ave,2)))

rh.ave <- rep(NA,4)
for(i in 1:4){
	rh.ave[i] <- mean(c(rh[i],rh[i+1]))
}
print(c("rh.ave"))
print(c(round(rh.ave,3)))

# Friction slope (energy loss gradient)
for(i in 2:5){
	sf[i - 1] <- (h.ead[i] - h.ead[i - 1])/ (xx[i] - xx[i - 1])
}
print(c("sf"))
print(c(round(sf,3)))

# WS slope
for(i in 2:5){
	sws[i - 1] <- (WSE[i] - WSE[i - 1])/ (xx[i] - xx[i - 1])
}
print(c("sws"))
print(c(round(sws,3)))

# Mannings n based on bed gradient
nn_so = rh.ave^(2/3) * so^0.5 / uu.ave
print(c("nn_so"))
print(c(round(nn_so,3)))

# Mannings n based on WSE
nn_sws = rh.ave^(2/3) * sws^0.5 / uu.ave
print(c("nn_sws"))
print(c(round(nn_sws,3)))

# Mannings n based on friction slope
nn_sf = rh.ave^(2/3) * sf^0.5 / uu.ave
print(c("nn_sf"))
print(c(round(nn_sf,3)))

f_sf <- 8 * 9.81 * rh.ave * sf / uu.ave^2
print(c("f_sf"))
print(c(round(f_sf,3)))

tau <- 1000*9.81*h.ave*sws
print(c("tau_sws"))
print(c(round(tau,0)))

omega <- 1000*9.81*Q*sws/wb.ave
print(c("omega_sws"))
print(c(round(omega,0)))


#q_temp <- rep(Q,5)

#o.ut <- cbind("probes"=probes,"x"=xx,"Q"=q_temp,"thead"=h.ead,"z"=z,"depths"=true_depths,"vel"=uu,"rh"=rh,"us_sf"=sf)
#o.ut

#write.csv(o.ut,"april_result_q_0.45.csv")

