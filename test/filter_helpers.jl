
# transition band is cosine. This allows satisfies the power complementary condition.
# cos( (yp-y)/(yp-ys) * π/2 ) for a falling edge. yp is pass band, ys is stop band.
function cosinetransitionlowpassfunction(y::T, yp::T, ys::T)::T where T
    abs_y = abs(y)

    if abs_y <= yp
        return one(T)
    elseif yp < abs_y <= ys
        #return cos( (yp-y)/(ys-yp)*π/2 )
        return cos( (yp-abs_y)/(ys-yp)*π/2 )
    end

    return zero(T)
end

function cosinetransitionhighpassfunction(x::T, xp::T, xs::T)::T where T
    abs_x = abs(x)

    if abs_x <= xs
        return zero(T)
    elseif xs < abs_x <= xp
        #return cos( (xp-x)/(xp-xs)*π/2 )
        return cos( (xp-abs_x)/(xp-xs)*π/2 )
    end

    return one(T)
end


# rs is rising edge's stopband.
# rp is rising edge's passband.
# fp is falling edge's passband.
# fs is falling edge's stopband.
function cosinetransitionbandpassfunction(y::T, rs, rp, fp, fs)::T where T
    abs_y = abs(y)

    # rising edge.
    if abs_y <= rs
        return zero(T)
    elseif rs < abs_y <= rp
        #return cos( (rp-y)/(rp-rs)*π/2 )
        return cos( (rp-abs_y)/(rp-rs)*π/2 )
    end

    # falling edge.
    if rp <= abs_y <= fp
        return one(T)
    elseif fp < abs_y <= fs
        #return cos( (fp-y)/(fs-fp)*π/2 )
        return cos( (fp-abs_y)/(fs-fp)*π/2 )
    end

    return zero(T)
end
