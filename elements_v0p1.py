# -*- coding: utf-8 -*-

########################################################################
#                                                                      #
#       Element classes for MAD-X-like tracking in lines with double   #
#       apertures.                                                     #
#       Author: L.S. Stoel                                             #
#       Version 0.1 - Work in progress                                 #
#                                                                      #
########################################################################

# Many simplifications assumed, e.g. linear optics and implicit
# assumption of dx/dt=px/ps=theta

# TODO Add more details on loss location

# TODO Implement reset() function for all elements, to return to initial
# parameters.
# TODO Implement noloss() function for all elements, to ignore aperture
# TODO Implement nofield() function for all elements, to set field to 0
# TODO Implement print() function for all elements, for debugging
# TODO Implement self.replacement for overcomplicated elements?
# TODO Add warnings for wrong arguments in constructor?
#      (eg negatice collimator thickness)

# TODO Add copy constuctor? Make constructor arguments optional?

import math

##################################################
#                                                #
#   Drift                                        #
#                                                #
##################################################

class Drift:
    """A 'drift' for accelerators."""
    def __init__(self, name, length, radius, offset_up=0, offset_down=0):
        self.name = name
        self.len = length
        self.r = radius
        self.offset_u = offset_up
        self.offset_d = offset_down

    def track(self, particle):
        # Does particle hit instantly?
        if self.r > 0:
            if (particle.x > self.offset_u+self.r
                    or particle.x < self.offset_u-self.r):
                particle.lost = self.name + '_start'
                return
        # Particle is within aperture!
        x_inc = self.len * particle.px
        # Does it hit downstream?
        if self.r > 0:
            if (particle.x + x_inc) > self.offset_d+self.r:
                s_hit_over_l = ((self.offset_u+self.r - particle.x)
                                / (x_inc + offset_u-offset_d))
                particle.s += self.len * s_hit_over_l
                particle.x += x_inc * s_hit_over_l
                particle.lost = self.name + '_down'
                return
            if (particle.x + x_inc) < self.offset_d-self.r:
                s_hit_over_l = ((self.offset_u-self.r - particle.x)
                                / (x_inc + offset_u-offset_d))
                particle.s += self.len * s_hit_over_l
                particle.x += x_inc * s_hit_over_l
                particle.lost = self.name + '_down'
                return
        # Particle made it out!
        particle.s += self.len
        particle.x += x_inc
        return

    def aperture(self, infty, s0):
        if self.r > 0:
            return [[[s0, self.r+self.offset_u],
                     [s0+self.len, self.r+self.offset_d],
                     [s0+self.len, self.r+self.offset_d+infty],
                     [s0, self.r+self.offset_u+infty]],
                    [[s0, self.offset_u-self.r],
                     [s0+self.len, self.offset_d-self.r],
                     [s0+self.len, self.offset_d-infty-self.r],
                     [s0, self.offset_u-infty-self.r]]]
        else:
            return []

##################################################
#                                                #
#   Kicker                                       #
#                                                #
##################################################

class Kicker:
    """A transverse kicker (dipole)."""
    def __init__(self, name, length, bendingangle, radius):
        self.name = name
        self.len = length
        self.an = bendingangle
        self.r = radius

    def track(self, particle):
        # Drift in disguise?
        if self.an == 0:
            Drift(self.name, self.len, self.r).track(particle)
            return
        # Actual kicker!
        if self.r > 0:
            # Does particle hit instantly?
            if particle.x > self.r or particle.x < (-1*self.r):
                particle.lost = self.name + '_start'
                return
            # Particle is within aperture!
            # Downstream hit on side the particle is bent away from?
            # (Solve an/len/2*s^2+px0*s+x0 == -sgn(an)*r)
            # (If solutions exist they are either both>0 or both<0)
            quadraticA = self.an/self.len/2
            quadraticB = particle.px
            quadraticC = particle.x + math.copysign(self.r, self.an)
            quadraticD = quadraticB**2 - 4*quadraticA*quadraticC
            if quadraticD > 0:
                hitdist = (-1*quadraticB - quadraticD**0.5) / (2*quadraticA)
                if hitdist > 0 and hitdist < self.len:
                    particle.s += hitdist
                    particle.x = -1*math.copysign(self.r, self.an)
                    particle.lost = self.name + '_down'
                    return
            # Downstream hit on side the particle is bent towards?
            # (Solve an/len/2*s^2+px0*s+x0 == sgn(an)*r)
            # (Guaranteed to have a solution >0 and one <0 because physics)
            quadraticC = particle.x - math.copysign(self.r, self.an)
            quadraticD = quadraticB**2 - 4*quadraticA*quadraticC
            hitdist = (-1*quadraticB + quadraticD**0.5) / (2*quadraticA)
            if hitdist < self.len:
                particle.s += hitdist
                particle.x = math.copysign(self.r, self.an)
                particle.lost = self.name + '_down'
                return
        # Particle made it out!
        particle.s += self.len
        particle.x += (particle.px + self.an/2) * self.len
        particle.px += self.an
        return

    def aperture(self, infty, s0):
        if self.r > 0:
            return [[[s0, self.r], [s0+self.len, self.r],
                     [s0+self.len, self.r+infty], [s0, self.r+infty]],
                    [[s0, 0-self.r], [s0+self.len, 0-self.r],
                     [s0+self.len, 0-infty-self.r], [s0, 0-infty-self.r]]]
        else:
            return []

##################################################
#                                                #
#   Quadrupole                                   #
#                                                #
##################################################

class Quadrupole:
    """A simple quadrupole, internal aperture check missing.

    Uses k=1/(Brho)*dBy/dx
    Aperture offsets define the central axis of the aperture, whereas
    field offset determines the neutral axis of the field. These are
    mostly meant to be used by the QuadHole class, but are also
    available to users.
    """

    def __init__(self, name, length, quad_k, radius, offset_field=0,
                 offset_aperture_up=0, offset_aperture_down=0):
        self.name = name
        self.len = length
        self.k = quad_k
        self.r = radius
        self.offset_au = offset_aperture_up
        self.offset_ad = offset_aperture_down
        self.offset_f = offset_field

    def track(self, particle):
        # Drift in disguise?
        if self.k == 0:
            Drift(self.name, self.len, self.r,
                  offset_up=self.offset_au,
                  offset_down=self.offset_ad).track(particle)
            return
        # Actual quadrupole!
        # Does particle hit instantly?
        if self.r > 0:
            if (particle.x > self.offset_au+self.r
                    or particle.x < self.offset_au-self.r):
                particle.lost = self.name + '_start'
                return
        # Particle is within aperture!
        xeff = particle.x - self.offset_f
        # Focussing quad?
        if self.k > 0:
            sk = self.k**0.5
            particle.s += self.len
            particle.x = (xeff*math.cos(sk*self.len)
                          + particle.px/sk*math.sin(sk*self.len)
                          + self.offset_f)
            particle.px = (-1.0*xeff*sk*math.sin(sk*self.len)
                           + particle.px*math.cos(sk*self.len))
        # Defocussing quad!
        else:
            sk = (-1.0*self.k)**0.5
            particle.s += self.len
            particle.x = (xeff*math.cosh(sk*self.len)
                          + particle.px/sk*math.sinh(sk*self.len)
                          + self.offset_f)
            particle.px = (xeff*sk*math.sinh(sk*self.len)
                           + particle.px*math.cosh(sk*self.len))
        # Downstream aperture check!
        if self.r > 0:
            if (particle.x > self.offset_ad+self.r
                    or particle.x < self.offset_ad-self.r):
                particle.lost = self.name + '_down'
                return

    def aperture(self, infty, s0):
        if self.r > 0:
            return [[[s0, self.r+self.offset_au], [s0+self.len, self.r+self.offset_ad],
                     [s0+self.len, self.r+self.offset_ad+infty], [s0, self.r+self.offset_au+infty]],
                    [[s0, self.offset_au-self.r], [s0+self.len, self.offset_ad-self.r],
                     [s0+self.len, self.offset_ad-infty-self.r], [s0, self.offset_au-infty-self.r]]]
        else:
            return []

##################################################
#                                                #
#   DoubleApDrift: drift with two apertures      #
#                                                #
##################################################

# TODO Make drift with angle to replace 0 thickness case?
# TODO Improve readability

class DoubleApDrift:
    """A drift with two separate apertures.

    Centered on the 'cirulating' aperture, with extraction aperture on
    the positive side.
    """
    def __init__(self, name, length, collpos_upstream, collpos_downstream,
                 coll_thickness, d_circulating, d_extraction):
        self.name = name
        self.len = length
        self.collpos_up = collpos_upstream
        self.collpos_down = collpos_downstream
        self.coll_thick = coll_thickness
        self.cdiam = d_circulating
        if d_circulating > 0:
            self.cdiam += coll_thickness/2
        self.ediam = d_extraction
        if d_extraction > 0:
            self.ediam += coll_thickness/2

    def track(self, particle):
        # Does particle hit instantly?
        # ...On the extraction side?
        if self.ediam > 0 and particle.x > (self.collpos_up + self.ediam):
            particle.lost = self.name + '_start_extr'
            return
        # ...On the collimator?
        if (particle.x < (self.collpos_up + self.coll_thick/2)
                and particle.x > (self.collpos_up - self.coll_thick/2)):
            particle.lost = self.name + '_start_coll'
            return
        # ...On the circulating side?
        if self.cdiam > 0 and particle.x < (self.collpos_up - self.cdiam):
            particle.lost = self.name + '_start_circ'
            return
        # Particle is within aperture!
        x_inc = self.len * particle.px
        # Circulating aperture?
        if particle.x < self.collpos_up:
            # Does it hit the downstream circulating aperture?
            if (self.cdiam > 0 and
                    (particle.x + x_inc) < (self.collpos_down - self.cdiam)):
                incfrac = ((particle.x - self.collpos_up + self.cdiam)
                           / (self.collpos_down - self.collpos_up - x_inc))
                particle.s += incfrac * self.len
                particle.x += incfrac * x_inc
                particle.lost = self.name + '_down_circ'
                return
            # Does it successfully exit from the circulating aperture?
            if ((particle.x + x_inc)
                    < (self.collpos_down - self.coll_thick/2)):
                particle.s += self.len
                particle.x += x_inc
                return
            # It reaches the collimator downstream!
            # Does it just go through an imaginary collimator?
            if not self.coll_thick > 0:
                # Does it hit the downstream extraction aperture?
                if (self.ediam > 0 and
                        (particle.x+x_inc) > (self.collpos_down + self.ediam)):
                    incfrac = ((particle.x - self.collpos_up - self.ediam)
                               / (self.collpos_down - self.collpos_up - x_inc))
                    particle.s += incfrac * self.len
                    particle.x += incfrac * x_inc
                    particle.lost = self.name + '_down_extr'
                    return
                # It successfully exits from the extraction aperture!
                particle.s += self.len
                particle.x += x_inc
                return
            # It is lost on the collimator!
            incfrac = ((particle.x - self.collpos_up)
                       / (self.collpos_down - self.collpos_up - x_inc))
            particle.s += incfrac * self.len
            particle.x += incfrac * x_inc
            particle.lost = self.name + '_down_coll_circ'
            return
        # Extraction aperture!
        # Does it hit the downstream extraction aperture?
        if (self.ediam > 0 and
                (particle.x + x_inc) > (self.collpos_down + self.ediam)):
            incfrac = ((particle.x - self.collpos_up - self.ediam)
                       / (self.collpos_down - self.collpos_up - x_inc))
            particle.s += incfrac * self.len
            particle.x += incfrac * x_inc
            particle.lost = self.name + '_down_extr'
            return
        # Does it successfully exit from the extraction aperture?
        if ((particle.x + x_inc)
                > (self.collpos_down + self.coll_thick/2)):
            particle.s += self.len
            particle.x += x_inc
            return
            # It reaches the collimator downstream!
            # Does it just go through an imaginary collimator?
            if not self.coll_thick > 0:
                # Does it hit the downstream circulating aperture?
                if (self.cdiam > 0 and
                        (particle.x+x_inc) < (self.collpos_down - self.cdiam)):
                    incfrac = ((particle.x - self.collpos_up + self.cdiam)
                               / (self.collpos_down - self.collpos_up - x_inc))
                    particle.s += incfrac * self.len
                    particle.x += incfrac * x_inc
                    particle.lost = self.name + '_down_circ'
                    return
                # It successfully exits from the circulating aperture!
                particle.s += self.len
                particle.x += x_inc
                return
            # It is lost on the collimator!
            incfrac = ((particle.x - self.collpos_up)
                       / (self.collpos_down - self.collpos_up - x_inc))
            particle.s += incfrac * self.len
            particle.x += incfrac * x_inc
            particle.lost = self.name + '_down_coll_extr'
            return

    def aperture(self, infty, s0):
        ans = [[[s0, self.collpos_up-self.coll_thick/2],
                 [s0+self.len, self.collpos_down-self.coll_thick/2],
                 [s0+self.len, self.collpos_down+self.coll_thick/2],
                 [s0, self.collpos_up+self.coll_thick/2]]]
        if self.ediam > 0:
            ans.extend([[[s0, self.collpos_up+self.ediam],
                 [s0+self.len, self.collpos_down+self.ediam],
                 [s0+self.len, self.collpos_down+self.ediam+infty],
                 [s0, self.collpos_up+self.ediam+infty]]])
        if self.cdiam > 0:
            ans.extend([[[s0, self.collpos_up-self.cdiam],
                 [s0+self.len, self.collpos_down-self.cdiam],
                 [s0+self.len, self.collpos_down-self.cdiam-infty],
                 [s0, self.collpos_up-self.cdiam-infty]]])
        return ans

##################################################
#                                                #
#   Septum                                       #
#                                                #
##################################################

# TODO have field region on either side possible

class Septum:
    """A septum with dipole extraction on positive side."""
    def __init__(self, name, length, bendingangle, bladepos_upstream,
                 bladepos_downstream, blade_thickness, d_circulating,
                 d_extraction):
        self.name = name
        self.len = length
        self.an = bendingangle
        self.bladepos_up = bladepos_upstream
        self.bladepos_down = bladepos_downstream
        self.blade_thick = blade_thickness
        self.cdiam = d_circulating
        if d_circulating > 0:
            self.cdiam += blade_thickness/2
        self.ediam = d_extraction
        if d_extraction > 0:
            self.ediam += blade_thickness/2

# TODO rewrite with track_circ and track_extr?
    def track(self, particle):
        # Double drift in disguise?
        if self.an == 0:
            temp_cdiam = 0
            temp_ediam = 0
            if self.cdiam > 0:
                temp_cdiam = self.cdiam-self.blade_thick/2
            if self.ediam > 0:
                temp_ediam = self.ediam-self.blade_thick/2
            DoubleApDrift(self.name, self.len, self.bladepos_up,
                          self.bladepos_down, self.blade_thick,
                          temp_cdiam, temp_ediam).track(particle)
            return
        # Actual septum!
        # Does particle hit instantly?
        # ...On the extraction side?
        if self.ediam > 0 and particle.x > (self.bladepos_up + self.ediam):
            particle.lost = self.name + '_start_extr'
            return
        # ...On the blade?
        if (particle.x < (self.bladepos_up + self.blade_thick/2)
                and particle.x > (self.bladepos_up - self.blade_thick/2)):
            particle.lost = self.name + '_start_blade'
            return
        # ...On the circulating side?
        if self.cdiam > 0 and particle.x < (self.bladepos_up - self.cdiam):
            particle.lost = self.name + '_start_circ'
            return
        # Particle is within aperture!
        # Circulating aperture?
        if particle.x < self.bladepos_up:
            x_inc = self.len * particle.px
            # Does it hit the downstream circulating aperture?
            if (self.cdiam > 0 and
                    (particle.x + x_inc) < (self.bladepos_down - self.cdiam)):
                incfrac = ((particle.x - self.bladepos_up + self.cdiam)
                           / (self.bladepos_down - self.bladepos_up - x_inc))
                particle.s += incfrac * self.len
                particle.x += incfrac * x_inc
                particle.lost = self.name + '_down_circ'
                return
            # Does it successfully exit from circulating aperture?
            if ((particle.x + x_inc)
                    < (self.bladepos_down - self.blade_thick/2)):
                particle.s += self.len
                particle.x += x_inc
                return
            # It reaches the blade downstream!
            incfrac = ((particle.x - self.bladepos_up - self.blade_thick/2)
                       / (self.bladepos_down - self.bladepos_up - x_inc))
            particle.s += incfrac * self.len
            particle.x += incfrac * x_inc
            # ... And hits it?
            if self.blade_thick > 0:
                particle.lost = self.name + '_down_blade_circ'
                return
            # ... And goes through the virtual blade!
            particle.update_history()
            templ = self.len - incfrac * self.len
            tempan = self.an - incfrac * self.an
            bladeposmid = particle.x
            # ... ... And hits the extraction aperture?
            if self.ediam > 0:
                quadraticA = tempan/templ/2
                quadraticB = (particle.px
                              - (self.bladepos_down - bladeposmid) / templ)
                quadraticC = particle.x - bladeposmid - self.ediam
                quadraticD = quadraticB**2 - 4*quadraticA*quadraticC
                hitdist = (-1*quadraticB + quadraticD**0.5) / (2*quadraticA)
                if hitdist < templ:
                    particle.s += hitdist
                    particle.x += (quadraticA * hitdist**2
                                   + particle.px * hitdist)
                    particle.px += tempan * hitdist / templ
                    particle.lost = self.name + '_down_extr'
                    return
            # ... ... And exits from the extraction aperture!
            particle.s += templ
            particle.x += (tempan/2 + particle.px) * templ
            particle.px += tempan
            return
        # Extraction aperture!
        quadraticA = self.an/self.len/2
        quadraticB = (particle.px
                      - (self.bladepos_down - self.bladepos_up) / self.len)
        quadraticC = particle.x - self.bladepos_up - self.blade_thick/2
        quadraticD = quadraticB**2 - 4*quadraticA*quadraticC
        # Does it reach the downstream blade?
        hitdist = (-1*quadraticB - quadraticD**0.5) / (2*quadraticA)
        if quadraticD > 0:
            if hitdist > 0 and hitdist < self.len:
                particle.s += hitdist
                particle.x += quadraticA * hitdist**2 + particle.px * hitdist
                particle.px += self.an * hitdist / self.len
                # ... And hit it?
                if self.blade_thick > 0:
                    particle.lost = self.name + '_down_blade_extr'
                    return
                # ... It goes through the virtual blade!
                particle.update_history()
                templ = self.len - hitdist
                x_inc = templ * particle.px
                bladeposmid = particle.x
                # ... ... And hits the circulating aperture?
                if (self.cdiam > 0 and
                        (particle.x+x_inc) < (self.bladepos_down-self.cdiam)):
                    incfrac = ((particle.x - bladeposmid + self.cdiam)
                               / (self.bladepos_down - bladeposmid - x_inc))
                    particle.s += incfrac * templ
                    particle.x += incfrac * x_inc
                    particle.lost = self.name + '_down_circ'
                    return
                # ... ... And exits from the extraction aperture!
                particle.s += templ
                particle.x += x_inc
                return
        # Does it hit the downstream extraction aperture?
        if self.ediam > 0:
            quadraticC = particle.x - self.bladepos_up - self.ediam
            quadraticD = quadraticB**2 - 4*quadraticA*quadraticC
            hitdist = (-1*quadraticB + quadraticD**0.5) / (2*quadraticA)
            if hitdist < self.len:
                particle.s += hitdist
                particle.x += quadraticA * hitdist**2 + particle.px * hitdist
                particle.px += self.an * hitdist / self.len
                particle.lost = self.name + '_down_extr'
                return
        # It successfully exits the extraction aperture!
        particle.s += self.len
        particle.x += (self.an/2 + particle.px) * self.len
        particle.px += self.an
        return

    def aperture(self, infty, s0):
        ans = [[[s0, self.bladepos_up-self.blade_thick/2],
                 [s0+self.len, self.bladepos_down-self.blade_thick/2],
                 [s0+self.len, self.bladepos_down+self.blade_thick/2],
                 [s0, self.bladepos_up+self.blade_thick/2]]]
        if self.ediam > 0:
            ans.extend([[[s0, self.bladepos_up+self.ediam],
                 [s0+self.len, self.bladepos_down+self.ediam],
                 [s0+self.len, self.bladepos_down+self.ediam+infty],
                 [s0, self.bladepos_up+self.ediam+infty]]])
        if self.cdiam > 0:
            ans.extend([[[s0, self.bladepos_up-self.cdiam],
                 [s0+self.len, self.bladepos_down-self.cdiam],
                 [s0+self.len, self.bladepos_down-self.cdiam-infty],
                 [s0, self.bladepos_up-self.cdiam-infty]]])
        return ans

##################################################
#                                                #
#   QuadHole: quadrupole with hole in yoke       #
#                                                #
##################################################

class QuadHole:
    """A quadrupole with a hole in the yoke.

    As in quadrupole, internal aperture check is missing and
    we use k=1/(Brho)*dBy/dx
    """

    def __init__(self, name, length, quad_k, hole_k, quad_radius,
                 hole_radius, hole_field_axis, hole_ap_axis_up,
                 hole_ap_axis_down):
        self.name = name
        self.len = length
        self.qk = quad_k
        self.hk = hole_k
        self.qr = quad_radius
        self.hr = hole_radius
        self.hfax = hole_field_axis
        self.haaxu = hole_ap_axis_up
        self.haaxd = hole_ap_axis_down

    def track(self, particle):
        # Is it in the hole?
        if particle.x > self.haaxu-self.hr and particle.x < self.haaxu+self.hr:
            Quadrupole(self.name+'_hole', self.len, self.hk, self.hr,
                       offset_field=self.hfax, offset_aperture_up=self.haaxu,
                       offset_aperture_down=self.haaxd).track(particle)
            return

        # Otherwise treat it like a normal quad.
        Quadrupole(self.name+'.circ', self.len, self.qk,
                   self.qr).track(particle)
        return

    def aperture(self, infty, s0):
        if self.qr <= 0 or self.hr <= 0:
            print("Aperture for QuadHole object ", self.name, " cannot be reliably drawn. Aperture omitted.")
            return []
        else:
            return [[[s0, 0-self.qr], [s0+self.len, 0-self.qr],
                     [s0+self.len, 0-self.qr-infty], [s0, 0-self.qr-infty]],
                    [[s0, self.qr], [s0+self.len, self.qr],
                     [s0+self.len, self.haaxd-self.hr], [s0, self.haaxu-self.hr]],
                    [[s0, self.haaxu+self.hr], [s0+self.len, self.haaxd+self.hr],
                     [s0+self.len, self.haaxd+self.hr+infty], [s0, self.haaxu+self.hr+infty]]]
