module datetime_utils

using Printf: @sprintf
using Dates: DateTime, Second, now, year
using ERFA: jd2cal


"""
Return how many seconds have passed since the
start of the current year
"""
function seconds_of_year()::Int64
    cur_time = now()
    start_of_year = DateTime(year(cur_time), 1, 1)

    since_start_of_year = cur_time - start_of_year
    # get rid of millisecond part
    since_start_of_year -= since_start_of_year % 1000

    return convert(Second, since_start_of_year).value
end



function caldate(jday::Float64)
    #==================================================

       This routine takes the modified Julian date and
       converts it to a date and time string.

       On Input:

          JDAY     modified Julian day (integer)

       On Output:

          DCHAR    date string (character)
          IDAY     day of the month (integer)
          IYDAY    year day (integer)
          IYEAR    year (integer)
          MONTH    month of the year (integer)
          TCHAR    time string (character)

       Calls:  GREGORIAN

    ==================================================#

    #--------------------------------------------------
    #  Define local data.
    #--------------------------------------------------
    #
    MCHAR = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    #            original from HOPS model

    OFFSET = 2440000

    #==================================================
       Begin executable code.
     =================================================#
    #
    #  Add offset to get true Julian date.

    julian = offset + trunc(Int64, jday)

    #fjulian = jday - int(jday) + 0.5
    fjulian = abs(jday - trunc(Int64, jday))

    if fjulian > 1.0
        julian += 1
        fjulian -= 1.0
    end

    #  Compute Gregorian date.
    gregdate = gregorian(julian)
    iday = gregdate.day
    month = gregday.month
    iyear = gregday.year
    iyday = gregday.yday

    #--------------------------------------------------
    #  Form date and time strings.
    #--------------------------------------------------

    hour = fjulian * 24.0
    ihour = int(hour)
    min = (hour - float(ihour)) * 60.0
    imin = int(min)
    isec = (min-imin) * 60.0

    dchar = @sprintf("%s%i%i", mchar[month], iday, iyear)
    tchar = @sprintf("%02i:%02i:%02i", ihour, imin, isec)
    return (
            dchar = dchar,
            iday  = iday,
            iyday = iyday,
            iyear = iyear,
            month = month,
            tchar = tchar,
           )
end

function julday(mm, id, iyyy)
    IGREG = 15 + 31 * (10 + 12 * 1582)
    iyyy == 0 && throw(DomainError("there is no year zero."))

    if iyyy < 0
        iyyy += 1
    end
    if mm > 2
        jy = iyyy
        jm = mm + 1
    else
        jy = iyyy - 1
        jm = mm + 13
    end
    julday = trunc(Int64, 365.25*jy) + trunc(Int64, 30.6001*jm) + id + 1720995

    if id + 31*(mm+12*iyyy) >= igreg
        ja = trunc(Int64, 0.01*jy)
        julday = julday + 2 - ja + int(0.25*ja)
    end

    return julday
end

function modjulianday(
        year::Int64, month::Int64, day::Int64, fracday::Float64
    )::Float64

    #   calculate the julian day from day, month, year and fraction of a day

    OFFSET = 2440000
    #                    original from HOPS model

    return julday(month, day, year) - OFFSET + fracday
end

struct GregorianDate
    day::Int32
    month::Int32
    year::Int32
    yday::Int32
end

const IYD = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
const IYDL = [1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367]
const IGREG = 2299161

# the following parameters are from wikipedia, which is from
# Richards 2013: Explanatory Supplement to the Astronomical Almanac, 3rd ed
# and Richards 1998: Mapping Time: The Calendar and its History
const YGREG::Int32 = 4716
const JGREG::Int32 = 1401
const MGREG::Int32 = 2
const NGREG::Int32 = 12
const RGREG::Int32 = 4
const PGREG::Int32 = 1461
const VGREG::Int32 = 3
const UGREG::Int32 = 5
const SGREG::Int32 = 153
const WGREG::Int32 = 2
const BGREG::Int32 = 274277
const CGREG::Int32 = -38
function gregorian(julian::Integer)::GregorianDate
    # all divisions here are meant to be integer divisions
    f = julian + JGREG + (((4 * julian + BGREG) ÷ 146097) * 3) ÷ 4 + CGREG
    e = RGREG * f + VGREG
    g = mod(e, PGREG) / RGREG
    h = UGREG * g + WGREG

    iday = mod(h, SGREG) ÷ UGREG + 1
    month = mod(h ÷ SGREG + MGREG, NGREG) + 1
    iyear = (e ÷ PGREG) - YGREG + (NGREG + MGREG - month) ÷ NGREG
end

end
