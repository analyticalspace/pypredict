import os
import time
from copy import copy
from cpredict import quick_find, quick_predict

try:
    basestring
except:
    basestring = str

def host_qth(path="~/.predict/predict.qth"):
    path = os.path.abspath(os.path.expanduser(path))
    try:
        with open(path) as qthfile:
            raw = [l.strip() for l in qthfile.readlines()]
            assert len(raw)== 4, "must match:\nname\nlat(N)\nlong(W)\nalt"%path
            return massage_qth(raw[1:])
    except Exception as e:
        raise RuntimeError("Unable to process qth '%s' (%s)"%(path, e))

def massage_tle(tle):
    try:
        # TLE may or may not have been split into lines already
        if isinstance(tle, basestring):
            tle = tle.rstrip().split('\n')
        assert len(tle) == 3, "TLE must be 3 lines, not %d: %s" % (len(tle), tle)
        return tle
        #TODO: print a warning if TLE is 'too' old
    except Exception as e:
        raise RuntimeError(e)

def massage_qth(qth):
    try:
        assert len(qth) == 3, "%s must consist of exactly three elements: (lat(N), long(W), alt(m))" % qth
        return (float(qth[0]), float(qth[1]), int(qth[2]))
    except ValueError as e:
        raise RuntimeError("Unable to convert '%s' (%s)" % (qth, e))
    except Exception as e:
        raise RuntimeError(e)

def observe(tle, qth, at=None):
    tle = massage_tle(tle)
    qth = massage_qth(qth)
    if at is None:
        at = time.time()
    return quick_find(tle, at, qth)

def transits(tle, qth, ending_after=None, ending_before=None):
    tle = massage_tle(tle)
    qth = massage_qth(qth)
    if ending_after is None:
        ending_after = time.time()
    ts = ending_after
    while True:
        transit = quick_predict(tle, ts, qth)
        t = Transit(tle, qth, start=transit[0]['epoch'], end=transit[-1]['epoch'])
        if (ending_before != None and t.end > ending_before):
            break
        if (t.end > ending_after):
            yield t
        # Need to advance time cursor so predict doesn't yield same pass
        ts = t.end + 60     #seconds seems to be sufficient

# Transit is a convenience class representing a pass of a satellite over a groundstation.
class Transit():
    def __init__(self, tle, qth, start, end):
        self.tle = massage_tle(tle)
        self.qth = massage_qth(qth)
        self.start = start
        self.end = end

        self.azimuth_start = observe(self.tle, self.qth, self.start)['azimuth']
        self.azimuth_end = observe(self.tle, self.qth, self.end)['azimuth']
        self.heading = self.find_heading()

        # Orbital velocity from PyPredict is in km/h
        self.peak_angular_rate = math.degrees(self.peak()['orbital_velocity'] / (self.peak()['slant_range']*3600))

    def find_heading(self):
        # First, convert compass to polar
        c2p = lambda x: 90-x if 90-x >= 0 else (360+(90-x)) % 360

        # Then, convert to radians
        polar_1 = math.radians(c2p(self.azimuth_start))
        polar_2 = math.radians(c2p(self.azimuth_end))

        # Then, convert polar to cartesian
        p2x = lambda x: (math.cos(x), math.sin(x))
        (x1, y1) = p2x(polar_1)
        (x2, y2) = p2x(polar_2)

        (x, y) = (x2 - x1, y2-y1)

        # Then, convert cartesian back to polar
        phi = math.degrees(math.atan2(y, x))

        # Finally, convert polar back to compass
        p2c = lambda x: 90-x if 90-x >= 0 else (360+(90-x)) % 360
        compass = p2c(phi)
        return compass

    # return observation within epsilon seconds of maximum elevation
    # NOTE: Assumes elevation is strictly monotonic or concave over the [start,end] interval
    def peak(self, epsilon=0.1):
        ts =  (self.end + self.start)/2
        step = (self.end - self.start)
        while (step > epsilon):
            step /= 4
            # Ascend the gradient at this step size
            direction = None
            while True:
                mid   = observe(self.tle, self.qth, ts)['elevation']
                left  = observe(self.tle, self.qth, ts - step)['elevation']
                right = observe(self.tle, self.qth, ts + step)['elevation']
                # Break if we're at a peak
                if (left <= mid >= right):
                    break
                # Ascend up slope
                slope = -1 if (left > right) else 1
                # Break if we stepped over a peak (switched directions)
                if direction is None:
                    direction = slope
                if direction != slope:
                    break
                # Break if stepping would take us outside of transit
                next_ts = ts + (direction * step)
                if (next_ts < self.start) or (next_ts > self.end):
                    break
                # Step towards the peak
                ts = next_ts
        return self.at(ts)

    # Return portion of transit above a certain elevation
    def above(self, elevation):
        return self.prune(lambda ts: self.at(ts)['elevation'] >= elevation)

    # Return section of a transit where a pruning function is valid.
    # Currently used to set elevation threshold, unclear what other uses it might have.
    # fx must either return false everywhere or true for a contiguous period including the peak
    def prune(self, fx, epsilon=0.1):
        peak = self.peak()['epoch']
        if not fx(peak):
            start = peak
            end = peak
        else:
            if fx(self.start):
                start = self.start
            else:
                # Invariant is that fx(right) is True
                left, right = self.start, peak
                while ((right - left) > epsilon):
                    mid = (left + right)/2
                    if fx(mid):
                        right = mid
                    else:
                        left = mid
                start = right
            if fx(self.end):
                end = self.end
            else:
                # Invariant is that fx(left) is True
                left, right = peak, self.end
                while ((right - left) > epsilon):
                    mid = (left + right)/2
                    if fx(mid):
                        left = mid
                    else:
                        right = mid
                end = left
        # Use copy to allow subclassing of Transit object
        pruned = copy(self)
        pruned.start = start
        pruned.end = end
        return pruned

    def duration(self):
        return self.end - self.start

    def at(self, t):
        if t < self.start or t > self.end:
            raise RuntimeError("time %f outside transit [%f, %f]" % (t, self.start, self.end))
        return observe(self.tle, self.qth, t)
