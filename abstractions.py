from typing import Dict, Set

import networkx as nx

class Tool:
    def __init__(self, name: str, version: str = None, container_uri: str = None) -> None:
        self.name: str = name
        self.version: str = version
        self.container_uri: str = container_uri

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
        return f"Dataflow(source={self.source.name}, destionation={self.destination.name}, data={[datum for datum in self.data]})"

class Workflow: 
    def __init__(self, name: str, steps: Set[Step], dataflows: Dataflow) -> None:
        self.name: str = name
        self.steps: Set[Step] = steps
        self.dataflows: Set[Dataflow] = dataflows

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

    def get_flattened_structure(self) -> nx.DiGraph:
        """ Get DAG representing flattened (i.e. no Subworkflow) structure of workflow. Mathematically equivalent to self.get_structure()

        Returns:
            nx.DiGraph: structure DAG.
        """
        g_structure: nx.DiGraph = nx.DiGraph()
        for connection in self.dataflows:
            source: Step = connection.source
            destination: Step = connection.destination
            g_structure.add_nodes_from([source, destination])
            g_structure.add_edge(source, destination)