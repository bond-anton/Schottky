class Dislocation(object):
    '''
    Linear defect in crystal
    '''

    def __init__(self, b, title):
        self.b = b
        self.title = title
        self.traps = []

    def add_trap(self, trap, concentration):
        self.traps.append([trap, concentration])

    def __str__(self, *args, **kwargs):
        description = 'Dislocation: ' + self.title + '\n'
        description += ' b: ' + str(self.b * 1e10) + ' Angstrom\n'
        for trap in self.traps:
            description += str(trap[0]) + '\n'
            description += ' Nl: ' + str(trap[1] * 1e-2) + ' 1/cm\n'
        return description
