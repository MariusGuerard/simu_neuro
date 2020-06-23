### Definition of the class of receptors mainly use for ORs.

class Receptor:
    """
    Receptor is used when the flag Receptor_flag = True in set_parameter.py
    """


    def __init__(self, NRec = 20, Nodor = 10):

        # Contains all the ORs

        self.ORs = []
