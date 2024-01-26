from logger import logger
from typing import Dict, Set

import networkx as nx

class Tool:
    def __init__(self, name: str, version: str = None, container_uri: str = None) -> None:
        self.name: str = name
        self.version: str = version
        self.container_uri: str = container_uri

    def __repr__(self) -> str:
        return f"Tool(name={self.name})"

class Datum:
    def __init__(self, uuid: str) -> None:
        self.uuid: str = uuid
        self.size: int = None

    def __repr__(self) -> str:
        return f"Datum(uuid={self.uuid[:5]})"

class Step:
    def __init__(self, name: str, is_tool: bool, tool: Tool = None) -> None:
        self.name: str = name
        self.outputs: Dict[Datum, str] = dict()
        self.inputs: Dict[Datum, str] = dict()
        self.is_tool: bool = is_tool
        self.tool: Tool = tool

        assert(self.check_valid()), f"{self} is_tool: {self.is_tool} but tool: {self.tool}"

    def __repr__(self) -> str:
        return f"Step(name={self.name}, inputs={self.inputs}, outputs={self.outputs}, tool={self.tool})"

    def register_output(self, datum: Datum, name: str) -> None:
        self.outputs[datum] = name

    def register_input(self, datum: Datum, name: str) -> None:
        self.inputs[datum] = name

    def check_valid(self) -> bool:
        return (
            (self.is_tool and self.tool)
            or (not self.is_tool and not self.tool)
        )

class Dataflow:
    def __init__(self, source: Step, destination: Step, data: Datum | Set[Datum]) -> None:
        self.source: Step = source
        self.destination: Step = destination
        self.data: Set[Datum] = (data if isinstance(data, Set) else set([data]))

    def add_data(self, datum: Datum) -> None:
        self.data.add(datum)

    def get_name_in_source(self, datum: Datum) -> str:
        return self.source.outputs[datum]

    def get_name_in_destination(self, datum: Datum) -> str:
        return self.source.inputs[datum]

    def __repr__(self) -> str:
        return f"Dataflow(source={self.source.name}, destination={self.destination.name}, data={[datum for datum in self.data]})"

class Workflow(Step): 
    def __init__(self, name: str, dataflows: Set[Dataflow]) -> None:
        self.name: str = name
        self.dataflows: Set[Dataflow] = dataflows
        Step.__init__(self, name=self.name, is_tool=False)
        self.steps: Set[Step] = self.get_steps()
        self.structure: nx.DiGraph = self.get_structure()

    def get_steps(self) -> Set[Step]:
        steps: Set = set()
        for flow in self.dataflows:
            steps.add(flow.source)
            steps.add(flow.destination)

        return steps

    def get_structure(self) -> nx.DiGraph:
        """ Get DAG representing structure of workflow.

        Returns:
            nx.DiGraph: structure DAG.
        """
        g_structure: nx.DiGraph = nx.DiGraph()
        for connection in self.dataflows:
            source: Step = connection.source
            destination: Step = connection.destination
            g_structure.add_nodes_from([source, destination])
            g_structure.add_edge(source, destination)

        return g_structure

    def get_flattened_structure(self) -> nx.DiGraph:
        """ Get DAG representing flattened (i.e. no Subworkflow) structure of workflow. Mathematically equivalent to self.get_structure()

        Returns:
            nx.DiGraph: structure DAG.
        """
        g_flat_structure: nx.DiGraph = self.get_structure()
        logger.info(f"{g_flat_structure.nodes}")
        for step in g_flat_structure.nodes:
            if isinstance(step, Workflow):
                inflows: Dict[Dataflow] = self.get_input_dataflows(self.structure, step)
                outflows: Dict[Dataflow] = self.get_output_dataflows(self.structure, step)
                print(f"{step} is a Workflow and has inflows={inflows}, outflows={outflows}")
                self.replace_inflows_destinations(self.structure, inflows, g_flat_structure)
                self.replace_outflows_sources(self.structure, outflows, g_flat_structure)

    def get_input_dataflows(self, step: Step) -> Set[Dataflow]:
        """ Get dataflows which have step as destination (i.e. step's input).

        Args:
            step (Step): target step.

        Returns:
            Set[Dataflow]: Dataflows with step as destination.
        """
        inflows: Set[Dataflow] = set()
        for flow in self.dataflows:
            if step is flow.destination:
                inflows.add(flow)

        return inflows

    def get_output_dataflows(self, step: Step) -> Set[Dataflow]:
        """ Get dataflows which have step as source (i.e. step's output).

        Args:
            step (Step): target step.

        Returns:
            Set[Dataflow]: Dataflows with step as source.
        """
        outflows: Set[Dataflow] = set()
        for flow in self.dataflows:
            if step is flow.source:
                outflows.add(flow)

        return outflows

    def flatten_step(self, dag: nx.DiGraph, step: Step, flattened_step: nx.DiGraph) -> nx.DiGraph:
        """ Flatten dag by replacing step with flattened_step. 

        Args:
            dag (nx.DiGraph): original dag.
            step (Step): step to be replaced.
            flattened_step (nx.DiGraph): flattened structure of step. 

        Returns:
            nx.DiGraph: flattened dag.
        """
        pass        

    def get_input_dataflows(self, dag: nx.DiGraph, step: Step) -> Set[Dataflow]:
        pass

    def get_output_dataflows(self, dag: nx.DiGraph, step: Step) -> Set[Dataflow]:
        pass

    def replace_inflows_destinations(self, original_structure: nx.DiGraph, dataflows: Set[Dataflow], target_structure: nx.DiGraph) -> None:
        pass


    def replace_outflows_sources(self, original_structure: nx.DiGraph, dataflows: Set[Dataflow], target_structure: nx.DiGraph) -> None:
        pass
