from __future__ import print_function

class Particle:
    """Particle-like object, used for tracking"""
    def __init__(self, x, px):
        self.s = 0
        self.x = x
        self.px = px
        self.lost = 'CIRCULATING'
        self.start = [0, x, px]
        self.history = [self.start]

    def state(self):
        return [self.s, self.x, self.px]

    def update_history(self):
        self.history.append(self.state())

    def print(self):
        print(" Particle state: " + str(self.state()) +
              "\n Lost?: " + self.lost +
              "\n Particle history: " + str(self.history) + "\n")
