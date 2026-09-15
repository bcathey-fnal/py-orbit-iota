"""
Module. The tune buffer node: hands the bunch to a monitor.tunebuffer
instance, which counts every particle's phase advance since the previous
node and, at the node marked as the first of the ring, closes the count of
the turn and records the bunch into its ring buffer.
"""
# import the base monitoring AccNode class
from ext.monitoring.monitoringAccNodes import Monitoring_Base_AccNode


class TuneBuffer_AccNode(Monitoring_Base_AccNode):
    """
    The subclass of the Monitoring_Base_AccNode for the tune buffer.
    """
    def __init__(self, tunebuffer, name="no name", isfirst=False):
        """
        Constructor. Creates the tune buffer node element.

        Parameters:
        - tunebuffer: the monitor.tunebuffer instance every node shares.
        - name: the node's name.
        - isfirst: whether this is the node at the start of the ring.
        """
        Monitoring_Base_AccNode.__init__(self, tunebuffer, "TuneBuffer", name)
        self.isfirst = isfirst

    def setFirst(self, isfirst):
        """
        Flag whether this is the node at the start of the ring.
        """
        self.isfirst = isfirst

    def __str__(self):
        return '{}: TuneBuffer, l = 0, first = {}, enable = {}'.format(
                self.getName(), self.isfirst, self.enable)

    def track(self, paramsDict):
        """
        Track the bunch through the node.
        """
        if(self.switcher != True): return
        bunch = paramsDict["bunch"]
        self.monitoring_calc.trackBunch(bunch, 1 if self.isfirst else 0)
