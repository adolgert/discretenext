"""
This runs a simulation using two different sampling techniques.
The model is specified in discrete time.

Usage:
    sample.py [-h] [-v] [-q] priority <filename> <runs> <stream>
    sample.py [-h] [-v] [-q] shifted <filename> <runs> <stream>
"""
import copy
import logging
import math
import time
import collections
import docopt
import numpy as np
import scipy
import scipy.stats
import networkx as nx
import h5py
import randomstate
import chemical
import pairing_heap
from gracefulinterrupthandler import GracefulInterruptHandler


logger=logging.getLogger("sample")


class GeometricTransition(object):
    def __init__(self, p):
        self.p=p
        self.dist=scipy.stats.geom(p=self.p)

    def Sample(self, rng):
        return rng.geometric(p=self.p)

    def Shifted(self, now, rng):
        return rng.binomial(n=1, p=self.p)==1


class UniformTransition(object):
    def __init__(self, a, b):
        assert(type(a)==int)
        assert(type(b)==int)
        assert(b>a)
        self.lims=[a, b]
        self.a=a
        self.b=b
        self.p=1.0/(b-a)

    def Sample(self, rng):
        return rng.randint(*self.lims)

    def Shifted(self, now,  rng):
        if now<self.a: return False
        assert(now<self.b)
        if self.b-now==1:
            return True
        logger.debug("UniformTransition b {0} now {0}".format(self.b, now))
        return rng.binomial(n=1, p=1.0/(self.b-now))==1


class Transition(object):
    def __init__(self, distribution, priority):
        self.distribution=distribution
        self.te=None
        self.priority=priority
        self.heap_node=None

    def Enable(self, now):
        self.te=now

    def Disable(self, now):
        self.te=None

    def Sample(self, now, rng):
        """
        This sampling method gets called the moment the transition is enabled,
        so it marks the enabling time, too. It asks when, in the future,
        the transition will fire.
        """
        self.te=now
        return now+self.distribution.Sample(rng)

    def SampleShifted(self, now, rng):
        """
        This sampling asks, given the current time, does this transition
        fire or not? It's a different sampling technique. The enabling
        time, self.te, will already be set.
        """
        logger.debug("Transition now {0} te {1}".format(now, self.te))
        return self.distribution.Shifted(now-self.te, rng)

    def Clone(self):
        return Transition(self.distribution, self.priority)



def CreateSISModel(N, step_max, transitions):
    survival=chemical.DiscreteSurvival(step_max)
    process=chemical.Process(survival)

    G=nx.complete_graph(N)
    node_to_idx=dict()
    initial_marking=dict()
    for ind_idx, pnode in enumerate(G.nodes_iter()):
        node_to_idx[pnode]=ind_idx

        for disease_state in ["S", "I"]:
            place_name=(ind_idx, disease_state)
            process.AddPlace(place_name, disease_state, 0)
            initial_marking[place_name]=0
    initial_marking[(0, "I")]=1
    for s_idx in range(1, N):
        initial_marking[(s_idx, "S")]=1

    for recover_idx in range(N):
        process.AddTransition(("R", recover_idx), "R",
            transitions["R"].Clone(),
            [((recover_idx, "I"), -1), ((recover_idx, "S"), 1)], 0)
    for source_n in G.nodes_iter():
        source_idx=node_to_idx[source_n]
        for target_n in G.neighbors(source_n):
            target_idx=node_to_idx[target_n]
            if source_idx != target_idx:
                process.AddTransition(("I", source_idx, target_idx), "I",
                    transitions["I"].Clone(),
                    [((source_idx, "I"), -1), ((source_idx, "I"), 1),
                    ((target_idx, "S"), -1), ((target_idx, "I"), 1)], 0)
    return process, initial_marking


def SamplePriority(model, initial_marking, step_cnt, summary, rng):
    """
    This sampling method draws for a future transition time
    at the moment a transition is enabled. This is written like
    Gibson and Bruck's Next Reaction method, except that it completely
    disobeys statistics by failing to draw from geometric distributions
    and then use a random variable transformation.
    """
    now=0
    last_step=-1
    heap=pairing_heap.pairing_heap()
    model.Reset(initial_marking, now)
    for first_t_name, first_t in model.AllEnabled():
        firing_time=first_t.Sample(now, rng)
        first_t.heap_node=heap.insert(
            (firing_time, first_t.priority, first_t_name))

    while not heap.empty():
        now, priority, who=heap.extract()
        should_be, was_enabled=model.Enabled(who)
        assert(should_be==was_enabled)
        assert(was_enabled)

        if now>step_cnt:
            break

        logger.debug("SamplePriority {0} {1}".format(now, who))
        model.Fire(who, now)

        disable, enable=model.AffectedTransitions()
        for dname, dtransition in disable:
            heap.delete(dtransition.heap_node)
            dtransition.heap_node=None
            dtransition.te=None
            model.Disable(dname, now)
        for ename, etransition in enable:
            efiring_time=etransition.Sample(now, rng)
            etransition.heap_node=heap.insert(
                (efiring_time, etransition.priority, ename))
            model.Enable(ename, now)
        if now!=last_step:
            summary[model.SummaryCounts()["I"]]+=1
            last_step=now

    model.FinishTiming(now)
    return now


def SampleShifted(model, initial_marking, step_cnt, summary, rng):
    """
    Think of Gillespie's First Reaction method. At every step,
    sample every enabled transition to see whether it will fire.
    Yes, this is incredibly slow.
    """
    now=0
    model.Reset(initial_marking, now)
    for fname, ftransition in model.AllEnabled():
        ftransition.te=now
    now=1

    while now<step_cnt:
        prioritized=collections.defaultdict(list)
        for first_t_name, first_t in model.AllEnabled():
            prioritized[first_t.priority].append((first_t_name, first_t))

        if not prioritized:
            break

        for priority_key in sorted(prioritized.keys()):
            for name, transition in prioritized[priority_key]:
                should_be, was_enabled=model.Enabled(name)
                assert(should_be==was_enabled)
                # It's possible a transition was disabled by another
                # transition scheduled for the same time.
                logger.debug("SampleShifted now {0} name {1}".format(now, name))
                if should_be and transition.SampleShifted(now, rng):
                    transition.te=None
                    model.Fire(name, now)

                    # How a transition affected the state of the system
                    # is usually calculated after a full sweep through all
                    # transitions, under the assumption that there are
                    # few or no conflicts. This assumption greatly reduces
                    # the order of computation, but putting this calculation
                    # here is assured to be correct in all cases.
                    disable, enable=model.AffectedTransitions()
                    for dname, dtransition in disable:
                        dtransition.te=None
                        model.Disable(dname, now)
                    for ename, etransition in enable:
                        etransition.te=now
                        model.Enable(ename, now)
                else:
                    pass
        now+=1
        summary[model.SummaryCounts()["I"]]+=1


    model.FinishTiming(now)
    return now


def WriteFile(filename, model, duration, summary):
    with GracefulInterruptHandler() as handler:
        out_data=h5py.File(filename, "w")
        grp_name="run{0:0>4d}".format(0)
        grp=out_data.create_group(grp_name)
        model.survival.write_hdf(grp)
        grp.create_dataset("duration", data=duration)
        grp.create_dataset("summary", data=summary)
        out_data.close()

if __name__ == "__main__":
    arguments = docopt.docopt(__doc__, version="sample 1.0")
    if arguments["-v"]:
        logging.basicConfig(level=logging.DEBUG)
    elif arguments["-q"]:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.INFO)

    rng=randomstate.prng.pcg64.RandomState(3333333, int(arguments["<stream>"]))

    filename=arguments["<filename>"]
    instance_cnt=int(arguments["<runs>"])
    N=10
    dt=0.01
    step_cnt=int(100/dt)
    # step_max is longest time it will take to fire a transition.
    step_max=int(10/dt)
    logger.info("step count {0}".format(step_cnt))
    priority={ "I" : 0, "R" : 1 }

    # The specification of distributions is discrete, but I'd like them
    # to behave similarly, not necessarily the same, as dt changes.
    # So we specify times in floating point and convert to integers.
    beta=0.5
    a=.2
    b=1.5
    transitions={
        "I" : Transition(GeometricTransition(beta*dt/(1+beta*dt)), priority["I"]),
        "R" : Transition(UniformTransition(round(a/dt), round(b/dt)), priority["R"])
    }

    model, initial_marking=CreateSISModel(N, step_max, transitions)

    if arguments["priority"]:
        sampler=SamplePriority
    elif arguments["shifted"]:
        sampler=SampleShifted

    time_limit_secs=10*60
    walltime=time.time()
    summary=np.zeros((N+1,), dtype=np.int)
    duration=np.zeros((step_cnt,), dtype=np.int)
    run_idx=0
    for i in range(instance_cnt):
        steps=sampler(model, initial_marking, step_cnt, summary, rng)
        if steps<step_cnt:
            duration[steps]+=1
        if time.time()-walltime > time_limit_secs:
            logger.info("Writing {0} to {1}".format(i, filename))
            WriteFile(filename, model, duration, summary)
            walltime=time.time()
            run_idx+=1
    WriteFile(filename, model, duration, summary)

    logger.info("Density\n{0}".format(summary))
    print_cnt=20
    end=np.where(duration>0)[0][-1]
    row_cnt=math.ceil(end/print_cnt)
    logger.info("Duration {0} out of total {1}".format(
        np.sum(duration), instance_cnt))
    for pr_idx in range(print_cnt):
        when=dt*pr_idx*row_cnt
        row_end=min((pr_idx+1)*row_cnt, end)
        amount=np.sum(duration[pr_idx*row_cnt:row_end])
        logger.info("{0:>8.2f} {1}".format(when, amount))

    survival=model.survival
    for dist_kind in ["I", "R"]:
        fire_cnt=np.sum(survival.fire[dist_kind])
        disable_cnt=np.sum(survival.disable[dist_kind])
        beyond=survival.beyond[dist_kind]
        logger.info("{0}: F {1} D {2} beyond {3}".format(dist_kind,
            fire_cnt, disable_cnt, beyond))