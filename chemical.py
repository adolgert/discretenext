import logging
import sys
import collections
import itertools
import numpy as np
import networkx as nx

logger=logging.getLogger("chemical")


class ContinuousSurvival(object):
    def __init__(self):
        self.fire=collections.defaultdict(list)
        self.disable=collections.defaultdict(list)

    def Fire(self, name, interval):
        self.fire[name].append(interval)

    def Disable(self, name, interval):
        self.disable[name].append(interval)

    def write_hdf(self, group):
        for k, v in self.fire.items():
            a=np.array(v, dtype=np.double)
            a.sort()
            group.create_dataset(str(k)+"fire", data=a)
            if k not in self.disable:
                self.disable[k]=list()

        for k, v in self.disable.items():
            a=np.array(v, dtype=np.double)
            a.sort()
            group.create_dataset(str(k)+"disable", data=a)


class DiscreteSurvival(object):
    def __init__(self, max_step):
        # The keys to these dictionaries are groups of events.
        # So that's Infection and Recovery for SIR.
        self.fire=collections.defaultdict(lambda: np.zeros((max_step,),
                dtype=np.int))
        self.disable=collections.defaultdict(lambda: np.zeros((max_step,),
                dtype=np.int))
        self.max_step=max_step
        # Beyond counts all firings or disablings that are greater than
        # max_step because they should be included in the total count.
        self.beyond=collections.defaultdict(int)

    def Fire(self, name, interval):
        """
        The name is a keyword for a group of events, such as Infection
        or recovery. The interval is an integer time between enabling
        and firing.
        """
        logger.debug("DiscreteSurvival::Fire {0}".format(interval))
        if interval<self.max_step:
            self.fire[name][interval]+=1
        else:
            logger.warn("DiscreteSurvival beyond fire")
            self.beyond[name]+=1

    def Disable(self, name, interval):
        """
        The name is a keyword for a group of events, such as Infection
        or recovery. The interval is an integer time between enabling
        and disabling.
        """
        logger.debug("DiscreteSurvival::Disable {0}".format(interval))
        if interval<self.max_step:
            self.disable[name][interval]+=1
        else:
            logger.warn("DiscreteSurvival beyond disable")
            self.beyond[name]+=1

    def write_hdf(self, group):
        for k, v in self.beyond.items():
            group.attrs["survivalbeyond"+k]=v
        for k, v in self.fire.items():
            group.create_dataset(str(k)+"fire", data=v)

        for k, v in self.disable.items():
            group.create_dataset(str(k)+"disable", data=v)


class Process(object):
    """
    This tracks the statistical process, with measurements of
    when transitions are enabled or disabled.
    Assumes that each transition has no memory, so it is enabled
    or disabled, never re-enabled.
    """
    def __init__(self, survival):
        self.G=nx.Graph()
        self.transition=dict()
        self.modified_place=set()
        self.transitions=dict()
        self.keyed_marking=collections.defaultdict(lambda: 0)
        self.survival=survival

    def AddPlace(self, name, key, now, value=0):
        self.G.add_node(name, v=value, time=now, key=key)
        # access with self.G.node[name] to get {"v" : value}.

    def AddTransition(self, name, key, transition, stoichiometry, now):
        self.transitions[name]=key
        self.G.add_node(name, t=transition, enabled=False, time=now)
        # Instead of using a multigraph, we define edge weights
        # as a list of weights, so "take one, give one" becomes
        # [1, -1] or [-1, 1].
        pw=collections.defaultdict(list)
        for place, weight in stoichiometry:
            logger.debug("AddTransition place {0} weight {1}".format(
                    place, weight))
            pw[place].append(weight)
        for p, w in pw.items():
            logger.debug("AddTransition p {0} w {1}".format(p, w))
            assert(self.G.node[p]["v"] is not None)
            self.G.add_edge(name, p, w=w)
            # Access with self.G.edge[a][b]["w"] to get weight.
        if self.Enabled(name):
            self.Enable(name, now)

    def Transition(self, name):
        return self.G.node[name]["t"]

    def FinishTiming(self, now):
        for name, key in self.transitions.items():
            if self.G.node[name]["enabled"]:
                self.Disable(name, now)

    def Reset(self, place_values, now):
        """
        Set the marking to the given place_values, with time at now.
        If a place isn't specified, it's set to zero.
        All transitions that can be enabled are enabled.
        place_values is a dictionary from place key to an integer
        token count. now is a time, which is integer for discrete
        and floating point for continuous.
        """
        self.keyed_marking=collections.defaultdict(lambda: 0)
        for place in [m for m in self.G.nodes() if "v" in self.G.node[m]]:
            if place in place_values:
                self.G.node[place]["v"]=place_values[place]
                self.G.node[place]["time"]=now
                self.keyed_marking[self.G.node[place]["key"]]+=place_values[place]
            else:
                self.G.node[place]["v"]=0
                self.G.node[place]["time"]=now
        for t in self.transitions.keys():
            enabled, was_enabled=self.Enabled(t)
            self.G.node[t]["time"]=now
            if enabled and not was_enabled:
                self._Enable(t, now)
            elif was_enabled and not enabled:
                self._Disable(t, now)

    def AllEnabled(self):
        transitions=list()
        for k, v in self.transitions.items():
            enabled_now, marked_enabled=self.Enabled(k)
            if enabled_now:
                transitions.append((k, self.G.node[k]["t"]))
        return transitions


    def Enabled(self, name):
        """
        Returns both whether transition should be enabled by the marking
        and whether it is currently marked as enabled.
        """
        assert("t" in self.G.node[name])
        for t, p in self.G.edges(name):
            for weight in self.G.edge[t][p]["w"]:
                if not "v" in self.G.node[p]:
                    logger.error("Enabled: node {0} {1}".format(p, self.G.node[p]))
                if weight+self.G.node[p]["v"]<0:
                    return (False, self.G.node[name]["enabled"])
        return (True, self.G.node[name]["enabled"])

    def DependencyGraph(self):
        """
        Assumes that each place has at most one token and calculates
        which transitions are mutually-exclusive.
        """
        dependency=nx.Graph()
        for place in [m for m in self.G.nodes() if "v" in self.G.node[m]
                and self.G.node[m]["v"]>0]:
            takes_a_token=list()
            borrows_a_token=list()
            for p, transition in self.G.edges(place):
                if self.G.node[transition]["enabled"]:
                    dependency.add_node(transition)
                    if sum(self.G.edge[p][transition]["w"]) is 0:
                        borrows_a_token.append(transition)
                    else:
                        takes_a_token.append(transition)
            logger.debug("DependencyGraph: p {0} take {1} borrow {2}".format(
                place, takes_a_token, borrows_a_token))
            for a, b in itertools.combinations(takes_a_token, 2):
                dependency.add_edge(a, b)
            for aa, bb in itertools.product(takes_a_token, borrows_a_token):
                dependency.add_edge(aa, bb)
        logger.debug("DependencyGraph exit {0}".format(dependency.nodes()))
        return dependency

    def Enable(self, name, now):
        self._Enable(name, now)
        return self.G.node[name]["t"]

    def Disable(self, name, now):
        self.survival.Disable(self.transitions[name],
                now-self.G.node[name]["time"])
        self._Disable(name, now)
        return self.G.node[name]["t"]

    def Fire(self, name, now):
        logger.debug("Fire: {0} {1}, {2}".format(now, name, self.Enabled(name)))
        if not self.G.node[name]["enabled"]:
            logger.error("Fire: of {0} at {1} but not enabled".format(
                name, now))
        self.survival.Fire(self.transitions[name],
                now-self.G.node[name]["time"])
        for t, p in self.G.edges(name):
            weight=sum(self.G.edge[t][p]["w"])
            # This marks a place as modified even if the weight is 0,
            # because the token moved.
            self.G.node[p]["v"]+=weight
            self.G.node[p]["time"]=now
            self.modified_place.add(p)
            self.keyed_marking[self.G.node[p]["key"]]+=weight
        self._Disable(name, now)
        return self.modified_place

    def SummaryCounts(self):
        return self.keyed_marking

    def AffectedTransitions(self):
        """
        Firing creates a list of modified places. This enables/disables
        transitions depending on those modified places.
        """
        enable=list()
        disable=list()
        neighbor_set=set()
        for p in self.modified_place:
            for t in nx.neighbors(self.G, p):
                neighbor_set.add(t)
        for neighbor in neighbor_set:
            (should_be, currently)=self.Enabled(neighbor)
            if should_be and not currently:
                enable.append((neighbor, self.G.node[neighbor]["t"]))
            elif currently and not should_be:
                disable.append((neighbor, self.G.node[neighbor]["t"]))
        return disable, enable

    def write_hdf(self, group):
        self.survival.write_hdf(group)


    def _Enable(self, name, now):
        n=self.G.node[name]
        if n["enabled"]:
            logger.error("_Enable {0} {1}".format(name, n["enabled"]))
            assert(not n["enabled"])
        logger.debug("_Enable {0}".format(name))
        n["enabled"]=True
        n["time"]=now

    def _Disable(self, name, now):
        """
        This does not record disabling of a transition in the survival.
        """
        n=self.G.node[name]
        if not n["enabled"]:
            logger.error("_Enable {0} {1}".format(name, n["enabled"]))
            assert(n["enabled"])
        logger.debug("_Disable {0}".format(name))
        n["enabled"]=False
        n["time"]=now


def CreateSIR(process, N, infect_transition, recover_transition):
    return CreateGraphSIR(process, infect_transition,
            recover_transition, nx.complete_graph(N))


def CreateGraphSIR(process, G, infect_transition, recover_transition):
    """
    Make SIR on a graph. If the graph is complete, then it's well-mixed
    SIR.
    """
    initial_marking=dict()
    node_to_idx=dict()
    N=len(G)
    logger.debug("CreateGraphSIR: nodes {0}".format(G.nodes()))
    for ind_idx, pnode in enumerate(G.nodes_iter()):
        node_to_idx[pnode]=ind_idx
        assert(ind_idx<N)
        for disease_state in ["S", "I", "R"]:
            place_name=(ind_idx, disease_state)
            process.AddPlace(place_name, disease_state, 0)
            initial_marking[place_name]=0
    logger.debug("CreateGraphSIR: ind_idx {0} N {1}".format(ind_idx, N))
    assert(ind_idx==N-1)
    initial_marking[(0, "I")]=1
    for s_idx in range(1, N):
        initial_marking[(s_idx, "S")]=1
    for recover_idx in range(N):
        process.AddTransition(("R", recover_idx), "R", recover_transition,
            [((recover_idx, "I"), -1), ((recover_idx, "R"), 1)], 0)
    for source_n in G.nodes_iter():
        source_idx=node_to_idx[source_n]
        logger.debug(G.neighbors(source_n))
        for target_n in G.neighbors(source_n):
            target_idx=node_to_idx[target_n]
            if source_idx!=target_idx:
                process.AddTransition(("I", source_idx, target_idx), "I",
                    infect_transition,
                    [((source_idx, "I"), -1), ((source_idx, "I"), 1),
                    ((target_idx, "S"), -1), ((target_idx, "I"), 1)], 0)
    logger.debug("CreateGraphSIR N {0} edges {1}".format(len(G), G.size()))
    return initial_marking
