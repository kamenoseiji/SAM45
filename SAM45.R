# This file stores functions to process SAM45 spectral data
#
source('date.R')	# to handle date into MJD
#
#-------- readSAM45 : a function to read SAM45 logging data
#
readSAM45 <- function(fname){
	tmpdata <- read.table(fname)
	rec_len <- dim(tmpdata)[2]  # elements in a record (= # of spectral channels + 3)
	spec <- as.matrix(tmpdata[4:rec_len])   # Extract spectal data
	spec <- spec / mean(spec)               # Normalize spectral data
	epoch <- tmpdata[,1]		    		# Epoch in YYYYMMDDHHMMSS
	year <- epoch %/% 10000000000
	doy  <- md2doy( year, (epoch %% 10000000000) %/% 100000000, (epoch %% 100000000) %/% 1000000)
	mjd  <- doy2fmjd( year, doy, (epoch %% 1000000) %/% 10000, (epoch %% 10000) %/% 100, epoch %% 100)
	scan <- data.frame(epoch, mjd, as.vector(tmpdata[2]), as.vector(tmpdata[3]))
	names(scan) <- c("epoch", "mjd", "IF", "scan")
	return(list(spec = t(spec), scan = scan))
}

#-------- BP : a function to calculate ON - smoothed OFF spectrum w/ Bandpass Calib (with a break)
#
BPcal_Predict <- function(SBCdata, variableBP, array_name){
	OnOff_flag <- which((SBCdata$scan$IF == array_name) & ( (SBCdata$scan$scan == 'ON') | (SBCdata$scan$scan == 'OFF')))
	SpecOnoff <- SBCdata$spec[273:3824,OnOff_flag]
	power <- colSums(SpecOnoff) / dim(SpecOnoff)[1]
	medPower <- median(power)
	powerExcess <- power - medPower
	BP <- numeric(nrow(SpecOnoff))
	for(scan_index in 1:dim(SpecOnoff)[2]){
		BP <- medPower* variableBP$BP + powerExcess[scan_index]* variableBP$ExcessSpec
		SpecOnoff[,scan_index] <- SpecOnoff[,scan_index] / BP
	}
	return(list(spec=SpecOnoff, scan=data.frame(epoch=SBCdata$scan$epoch[OnOff_flag], IF=SBCdata$scan$IF[OnOff_flag], scan=SBCdata$scan$scan[OnOff_flag])))
}
		

#-------- SBC_Spec_Break : a function to calculate ON - smoothed OFF spectrum w/ Bandpass Calib (with a break)
#
SBC_Spec_Break <- function( SBCdata, nsw, array_name, breakCH=1776 ){
	# SBCdata : Input spectral data
	# nsw     : Smoothing Window (number of channels)
	#---- ON and OFF scans
	SBC_ON <- SBCdata$spec[,which(SBCdata$scan$IF == array_name & (SBCdata$scan$scan=='ON'))]
	SBC_OFF <- SBCdata$spec[,which(SBCdata$scan$IF == array_name & (SBCdata$scan$scan=='OFF'))]
	
	#---- Spiline Smoothing for OFF scans
	SBC_spline <- matrix(nrow=nrow(SBC_ON), ncol=ncol(SBC_OFF))
	for(index in 1:ncol(SBC_OFF)){
		SBC_spline[,index] <- c(predict(smooth.spline(SBC_OFF[1:breakCH, index], all.knots=TRUE, df=(breakCH %/% nsw)), 1:breakCH)$y, predict(smooth.spline(SBC_OFF[(breakCH+1):dim(SBC_OFF)[1], index], all.knots=TRUE, df=((dim(SBC_OFF)[1]-breakCH) %/% nsw)), 1:(dim(SBC_OFF)[1]-breakCH))$y)
	}
	return( (rowSums(SBC_ON)/dim(SBC_ON)[2]) / (rowSums(SBC_spline)/dim(SBC_OFF)[2]) - 1 )
}

#-------- SBC_spec : a function to calculate ON - smoothed OFF spectrum
#
SBC_spec <- function( SBCdata, nsw, array_name ){
	# SBCdata : Input spectral data
	# nsw     : Smoothing Window (number of channels)
	#---- ON and OFF scans
	SBC_ON <- SBCdata$spec[,which(SBCdata$scan$IF == array_name & (SBCdata$scan$scan=='ON'))]
	SBC_OFF <- SBCdata$spec[,which(SBCdata$scan$IF == array_name & (SBCdata$scan$scan=='OFF'))]
	
	#---- Spiline Smoothing for OFF scans
	SBC_spline <- matrix(nrow=nrow(SBC_ON), ncol=ncol(SBC_OFF))
	for(index in 1:ncol(SBC_OFF)){
		SBC_spline[,index] <- predict(smooth.spline(SBC_OFF[, index], all.knots=TRUE, df=(dim(SBC_OFF)[1] %/% nsw)), 1:dim(SBC_OFF)[1])$y
	}
	return( (rowSums(SBC_ON)/dim(SBC_ON)[2]) / (rowSums(SBC_spline)/dim(SBC_OFF)[2]) - 1 )
}

#-------- SBC_spec : a function to calculate ON - smoothed OFF spectrum
#
SBC_BPspec <- function( SBCdata, nsw, array_name ){
	# SBCdata : Input spectral data
	# nsw     : Smoothing Window (number of channels)
	#---- ON and OFF scans
	SBC_ON <- SBCdata$spec[,which(SBCdata$scan$IF == array_name & (SBCdata$scan$scan=='ON'))]
	SBC_OFF <- SBCdata$spec[,which(SBCdata$scan$IF == array_name & (SBCdata$scan$scan=='OFF'))]
	#---- Make BP
	BP_raw <- rowSums(SBC_OFF) / dim(SBC_OFF)[2]	# Integrate all OFF scans
	BP <- predict(smooth.spline(BP_raw, all.knots=TRUE, df=(dim(SBC_OFF)[1] %/% 4)), 1:dim(SBC_OFF)[1])$y
	#---- BP Calib
	SBC_ON <- SBC_ON / BP; SBC_OFF <- SBC_OFF / BP
	
	#---- Spiline Smoothing for OFF scans
	SBC_spline <- matrix(nrow=nrow(SBC_OFF), ncol=ncol(SBC_OFF))
	for(index in 1:ncol(SBC_OFF)){
		SBC_spline[,index] <- predict(smooth.spline(SBC_OFF[, index], all.knots=TRUE, df=(dim(SBC_OFF)[1] %/% nsw)), 1:dim(SBC_OFF)[1])$y
	}
	return( (rowSums(SBC_ON)/dim(SBC_ON)[2]) / (rowSums(SBC_spline)/dim(SBC_OFF)[2]) - 1 )
}

#-------- STD_spec : a function to calculate conventional ON - OFF
#
STD_spec <- function( STDdata, array_name ){
	STD_ON  <- STDdata$spec[,which(STDdata$scan$IF == array_name & (STDdata$scan$scan=='ON'))]
	STD_OFF  <- STDdata$spec[,which(STDdata$scan$IF == array_name & (STDdata$scan$scan=='OFF'))]
	return( (rowSums(STD_ON)/dim(STD_ON)[2]) / ( rowSums(STD_OFF)/dim(STD_OFF)[2]) - 1 )
}

#-------- SBC_Cal : a function to produece smoothed-bandpass-calibrated spectrum
SBC_Cal <- function( rawSpec, BP, time_integ, node_ch ){
	# rawSpec : raw spectrum before calibration / integration / smoothing
	# BP      : template bandpass
	# time_integ : integration period
	# node_ch : Smoothing range
	#
	numSpec <- dim(rawSpec)[1]; numTime <- dim(rawSpec)[2]    # 周波数点数と時間点数
	caledSpec <- t(apply(rawSpec/BP, 1, bunch, lag=time_integ)) / numSpec   # BP補正したスペクトルを時間積分。
	Power_time <- colSums(caledSpec) / numSpec                # 平均電力の時間変化
	caledSpec <- t(t(caledSpec) / Power_time)                 # 平均電力で正規化
	SplineBP <- matrix(nrow=nrow(caledSpec), ncol=ncol(caledSpec))  # 出力の行列を作成
	for(scan_index in 1:ncol(SplineBP)){
		SplineBP[,scan_index] <- predict(smooth.spline(caledSpec[2:numSpec, scan_index], all.knots=TRUE, df=(numSpec %/% node_ch)), 1:numSpec)$y # OFF点スペクトルを周波数方向に平滑化
	}
	return( SplineBP )
}

#-------- OnOffSpecResult : a function to produece ON-OFF Scan spectrum
OnoffSpecResult <- function(rawSpec, scanSec = 32, integScan = 16){
	# rawSpec : raw spectrum before calibration / integration / smoothing
	# scanSec : integration time of each scan
	# integScan : Number of scans to be integrated
	
	#-------- convensional ON-OFF scans
	numSpec <- dim(rawSpec)[1]; numTime <- dim(rawSpec)[2]	# 周波数点数と時間点数
	integSpec <- t(apply(rawSpec, 1, bunch, lag=scanSec))	# scanSec ごとにスペクトルを時間積分。
	numScan <- dim(integSpec)[2] %/% (2* integScan)			# Number of whole scans
	on_pattern <- rep(c(0,1), integScan); off_pattern <- rep(c(1,0), integScan)
	spec <- matrix( nrow=numSpec, ncol=numScan)			# スペクトル結果を格納するmatrix
	spec_bl <- matrix( nrow=numSpec, ncol=numScan)			# スペクトル結果を格納するmatrix
	spec_sd <- numeric(numScan); spec_bl_sd <- numeric(numScan)	# SDを格納するベクトル
	for( scan_index in 1:numScan ){
		integ_start <- 2* integScan* (scan_index - 1) + 1
		integ_end   <- 2* integScan* scan_index
		cat(sprintf(" sample %d : integ = %d - %d\n", scan_index, integ_start, integ_end))
		scanSpec <- integSpec[, integ_start:integ_end]
		on_spec  <- rowSums(t(t(scanSpec)* on_pattern)) / sum(on_pattern)  # ON点をスキャンパターンに合わせて積分
		off_spec <- rowSums(t(t(scanSpec)* off_pattern))/ sum(off_pattern) # OFF点をスキャンパターンに合わせて積分
		spec[,scan_index] <- (on_spec - off_spec)/off_spec
		spec_sd[scan_index] <- sd(spec[,scan_index])
		plot(spec[,scan_index], type='l')
		
		#-------- Baseline Subtraction
		spec_bl[,scan_index] <- spec[,scan_index] - predict(smooth.spline(spec[,scan_index], all.knots=TRUE, df=(numSpec %/% 64)), (1:numSpec))$y  # baseline差引き
		spec_bl_sd[scan_index] <- sd(spec_bl[,scan_index])
		plot(spec_bl[,scan_index], type='l')
	}
	return( list(spec = spec, bl_spec =  spec_bl, sd_spec = spec_sd, sd_bl_spec = spec_bl_sd))
}


SplineSpecResult <- function(rawSpec, BP, scanSec = 8, OnOffRatio = 7, integScan = 8, node_ch=128){
	# rawSpec : raw spectrum before calibration / integration / smoothing
	# BP      : Template Bandpass
	# scanSec : Off点1回の積分時間
	# OnOffRatio : On点時間とOff点時間の比
	# integScan : (On + Off)のセットを何回積分するか
	# node_ch : Spline平滑化の幅
	numSpec <- dim(rawSpec)[1]; numTime <- dim(rawSpec)[2]		# 周波数点数と時間点数
	numScan <- numTime %/% (scanSec * (OnOffRatio + 1) * integScan)	# Number of whole scans
	#-------- BPによる補正
	caledSpec <- t(apply(rawSpec, 1, bunch, lag=scanSec))/BP	# スペクトルを時間積分してBP補正
	Power_time <- colSums(caledSpec)							# 平均電力の時間変化
	caledSpec <- t(t(caledSpec) / Power_time)					# 平均電力で正規化
	
	#-------- 平滑化したスペクトル
	SplineBP <- matrix(nrow=nrow(caledSpec), ncol=ncol(caledSpec))  # 平滑化したスペクトル格納用
	spec_sd <- numeric(numScan); spec_bl_sd <- numeric(numScan)	# SDを格納するベクトル
	for(scan_index in 1:ncol(SplineBP)){
		SplineBP[,scan_index] <- predict(smooth.spline(caledSpec[2:numSpec, scan_index], all.knots=TRUE, df=(numSpec %/% node_ch)), 1:numSpec)$y # OFF点スペクトルを周波数方向に平滑化
	}
	
	#-------- スキャンパターンに合わせて積分
	off_pattern <- rep(c(1, rep(0,OnOffRatio)), integScan)
	on_pattern  <- rep(c(0, rep(1,OnOffRatio)), integScan)
	spec <- matrix( nrow=numSpec, ncol=numScan)
	spec_bl <- matrix( nrow=numSpec, ncol=numScan)
	for( scan_index in 1:numScan ){
		integ_start <- (OnOffRatio + 1)* integScan* (scan_index - 1) + 1
		integ_end   <- (OnOffRatio + 1)* integScan* scan_index
		cat(sprintf(" sample %d : integ = %d - %d\n", scan_index, integ_start, integ_end))		
		on_spec  <- rowSums(t(t(caledSpec[, integ_start:integ_end])* on_pattern)) / sum(on_pattern)  # ON点をスキャンパターンに合わせて積分
		off_spec <- rowSums(t(t(SplineBP[, integ_start:integ_end])* off_pattern))/ sum(off_pattern) # OFF点をスキャンパターンに合わせて積分
		spec[,scan_index] <- (on_spec - off_spec)/off_spec
		spec_sd[scan_index] <- sd(spec[,scan_index])
		plot(spec[,scan_index], type='l')
	
		#-------- Baseline Subtraction
		spec_bl[,scan_index] <- spec[,scan_index] - predict(smooth.spline(spec[,scan_index], all.knots=TRUE, df=(numSpec %/% 64)), (1:numSpec))$y  # baseline差引き
		spec_bl_sd[scan_index] <- sd(spec_bl[,scan_index])
		plot(spec_bl[,scan_index], type='l')
	}
	return( list(spec = spec, bl_spec =  spec_bl, sd_spec = spec_sd, sd_bl_spec = spec_bl_sd))
}

SplineVarSpecResult <- function(rawSpec, BP, varBP, medPower, excessPower, scanSec = 8, OnOffRatio = 7, integScan = 8, node_ch=128){
	# rawSpec : raw spectrum before calibration / integration / smoothing
	# BP      : Template Bandpass
	# varBP   : Bandpass variable component
	# medPower: Power as a function of time
	# excessPower : Excess component
	# scanSec : Off点1回の積分時間
	# OnOffRatio : On点時間とOff点時間の比
	# integScan : (On + Off)のセットを何回積分するか
	# node_ch : Spline平滑化の幅
	numSpec <- dim(rawSpec)[1]; numTime <- dim(rawSpec)[2]		# 周波数点数と時間点数
	numScan <- numTime %/% (scanSec * (OnOffRatio + 1) * integScan)	# Number of whole scans
	#-------- BPによる補正
	caledTemp <- matrix(nrow=numSpec, ncol=numTime)				# BP補正作業用のmatrix
	for(time_index in 1:numTime){								# variable BPの補正
		caledTemp[,time_index] <- rawSpec[,time_index] / (medPower[time_index]* BP + excessPower[time_index]* varBP)
	}
	caledSpec <- t(apply(caledTemp, 1, bunch, lag=scanSec))		# スペクトルを時間積分
	Power_time <- colSums(caledSpec)							# 平均電力の時間変化
	caledSpec <- t(t(caledSpec) / Power_time)					# 平均電力で正規化
	
	#-------- 平滑化したスペクトル
	SplineBP <- matrix(nrow=nrow(caledSpec), ncol=ncol(caledSpec))  # 平滑化したスペクトル格納用
	spec_sd <- numeric(numScan); spec_bl_sd <- numeric(numScan)	# SDを格納するベクトル
	for(scan_index in 1:ncol(SplineBP)){
		SplineBP[,scan_index] <- predict(smooth.spline(caledSpec[2:numSpec, scan_index], all.knots=TRUE, df=(numSpec %/% node_ch)), 1:numSpec)$y # OFF点スペクトルを周波数方向に平滑化
	}
	
	#-------- スキャンパターンに合わせて積分
	off_pattern <- rep(c(1, rep(0,OnOffRatio)), integScan)
	on_pattern  <- rep(c(0, rep(1,OnOffRatio)), integScan)
	spec <- matrix( nrow=numSpec, ncol=numScan)
	spec_bl <- matrix( nrow=numSpec, ncol=numScan)
	for( scan_index in 1:numScan ){
		integ_start <- (OnOffRatio + 1)* integScan* (scan_index - 1) + 1
		integ_end   <- (OnOffRatio + 1)* integScan* scan_index
		cat(sprintf(" sample %d : integ = %d - %d\n", scan_index, integ_start, integ_end))		
		on_spec  <- rowSums(t(t(caledSpec[, integ_start:integ_end])* on_pattern)) / sum(on_pattern)  # ON点をスキャンパターンに合わせて積分
		off_spec <- rowSums(t(t(SplineBP[, integ_start:integ_end])* off_pattern))/ sum(off_pattern) # OFF点をスキャンパターンに合わせて積分
		spec[,scan_index] <- (on_spec - off_spec)/off_spec
		spec_sd[scan_index] <- sd(spec[,scan_index])
		plot(spec[,scan_index], type='l')
	
		#-------- Baseline Subtraction
		spec_bl[,scan_index] <- spec[,scan_index] - predict(smooth.spline(spec[,scan_index], all.knots=TRUE, df=(numSpec %/% 64)), (1:numSpec))$y  # baseline差引き
		spec_bl_sd[scan_index] <- sd(spec_bl[,scan_index])
		plot(spec_bl[,scan_index], type='l')
	}
	return( list(spec = spec, bl_spec =  spec_bl, sd_spec = spec_sd, sd_bl_spec = spec_bl_sd))
}

#-------- Median Window Filter
mwf <- function( x, window ){
	# x      : Input vector of values
	# window : window width
	num_x <- length(x)	# Length of the input vector
	y <- numeric(num_x)	# Output smoothed vector
	leg <- (window - 1) %/% 2	# Half length of window
	if(leg <= 1){ return(x)}	# Too narrow window

	for( index in 1:leg ){	y[index] <- median(x[1 : (index + leg)]) }	# Head of the seaquence
	for( index in (leg + 1):(num_x - leg)){	y[index] <- median(x[(index - leg) : (index + leg)])}	# middle of the sequence
	for( index in (num_x - leg + 1):num_x){ y[index] <- median(x[(index - leg) : num_x])}			# tail of the sequence
	return(y)
}

#-------- Extract spectrum of variable component
burstSpec <- function( rawSpec ){
	numSpec <- dim(rawSpec)[1]; numTime <- dim(rawSpec)[2]		# 周波数点数と時間点数
	power <- colSums(rawSpec)			# Power variation
	medPower <- mwf( power, 301 )		# Median Window Filter w/ 301-sec range
	powerExcess <- power - medPower		# Power excess
	
	burstIndex <- burstIndex <- which( powerExcess > 0.01)	# Filter for burst state
	steadyIndex <- which( powerExcess < 0.01)				# Filter for steady state
	BP <- rowSums(rawSpec[,steadyIndex]) / length(steadyIndex)	# make BP within steady state
	ExcessSpec <- as.vector(rawSpec[,burstIndex] %*% powerExcess[burstIndex]) / sum(powerExcess[burstIndex])
	ExcessSpec <- ExcessSpec / BP - 1
	node_ch <- which.min(allanvar(ExcessSpec)[1:128])
	return(predict(smooth.spline(ExcessSpec, all.knots=TRUE, df=(numSpec %/% node_ch)), 1:numSpec)$y)
}

#-------- Solve variable spectra into static and variable components 
variableBP <- function( rawSpec ){
	numSpec <- dim(rawSpec)[1]; numTime <- dim(rawSpec)[2]		# 周波数点数と時間点数
	power <- colSums(rawSpec) / dim(rawSpec)[1]					# Power variation
	medPower <- mwf( power, 301 )		# Median Window Filter w/ 301-sec range
	powerExcess <- power - medPower		# Power excess
	
	BP <- numeric(numSpec); ExcessSpec <- numeric(numSpec)
	for(index_ch in 1:numSpec){
		fit <- lm(formula = rawSpec[index_ch,] ~ medPower + powerExcess + 0)
		BP[index_ch] <- fit$coefficients[1]
		ExcessSpec[index_ch] <- fit$coefficients[2]
#		ExcessSpec[index_ch] <- 0
	}
	
	return(list(power = power, medPower = medPower, powerExcess = powerExcess, BP = BP, ExcessSpec = ExcessSpec))
}
