MJD_1901 <- 15384		# MJD at Y1901
DAY_PER_YEAR <- 365
DAY_PER_4YEAR <- 1461
SEC_PER_DAY <- 86400
DOY_MON <- function(month) c(0, 31, 59, 90, 120, 151, 181, 212, 242, 273, 303, 334)[month]
DOY_MON_LEAP <- function(month) c(0, 31, 60, 91, 121, 152, 182, 213, 243, 274, 304, 335)[month]

#-------- Calculate Modified Julian Date from Year, Month, and Day
day2mjd <- function(year, month, day){
	doy <- md2doy(year, month, day)
	doy2mjd(year, doy)
}

#-------- Calculate Modified Julian Date from Year and Day of Year
doy2mjd <- function( year, doy ){
	if((year < 1901) || (year > 2099)){ return(-1)}
	year <- year - 1901
	mjd <- year %/%4 * DAY_PER_4YEAR + year%%4 * DAY_PER_YEAR + doy + MJD_1901
	return(mjd)
}

#-------- Calculate (fractional) Modified Julian Date
doy2fmjd <- function(year, doy, hour, min, sec){
	if((year < 1901) || (year > 2099)){ return(-1)}
	year <- year - 1901
	mjd <- mjd <- year %/%4 * DAY_PER_4YEAR + year%%4 * DAY_PER_YEAR + doy + MJD_1901
	sod <- (hour*60 + min)*60 + sec
	mjd <- mjd + sod/SEC_PER_DAY
	return(mjd)
}

#-------- Calculate Day of year from Month and date
md2doy <- function(year, month, date){
	is_leap <- ((year%%4 == 0) && (month > 2))	# Leap Year Frag
	DOY_MON(month) + date + is_leap
}
