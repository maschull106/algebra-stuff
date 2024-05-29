from __future__ import annotations
from .hilbert_scheme import *
from typing import Tuple


class YoungDiagramIdeals:
    def __init__(self, base: DoubleNestedHilbertScheme, nested_ideals: List[List[PolyRingIdeal]], R: PolyRing = None):
        self.base = base
        self.diagram = base.diagram
        self.nested_ideals = nested_ideals
        if R is None:
            R = infer_poly_ring()
        self.R = R
        self.sanity_check()
    
    def sanity_check(self):
        if any(self.diagram[i+1] > self.diagram[i] for i in range(len(self.diagram)-1)):
            raise ValueError
        if len(self.nested_ideals) != len(self.diagram):
            raise ValueError
        if any(len(row) != row_len for row, row_len in zip(self.nested_ideals, self.diagram)):
            return ValueError
        if any(ideal.base != self.R for row in self.nested_ideals for ideal in row):
            raise ValueError
        
        inclusion_constraint = lambda seq: all(seq[i].contains_ideal(seq[i+1]) for i in range(len(seq)-1))
        if not all(map(inclusion_constraint, self.rows())):
            raise ValueError
        if not all(map(inclusion_constraint, self.columns())):
            raise ValueError
    
    def __getitem__(self, pos: Tuple[int, int]):
        """row i, column j"""
        i, j = pos
        return self.nested_ideals[i][j]
    
    def __iter__(self):
        for row in self.nested_ideals:
            for ideal in row:
                yield ideal
    
    def size(self) -> int:
        return sum(map(len, self.nested_ideals))
    
    def row(self, i: int) -> List[PolyRingIdeal]:
        return self.nested_ideals[i]
    
    def column(self, i: int) -> List[PolyRingIdeal]:
        column = []
        row_ind = 0
        while row_ind < len(self.diagram) and self.diagram[row_ind] >= i+1:
            column.append(self.nested_ideals[row_ind][i])
            row_ind += 1
        return column
    
    def rows(self) -> List[List[PolyRingIdeal]]:
        return self.nested_ideals
    
    def columns(self) -> List[List[PolyRingIdeal]]:
        return [self.column(i) for i in range(self.diagram[0])]
    
    def index_mapping(self, i: int, j: int) -> int:
        """convert index (i, j) to the corresponding index if all the rows were concatenated"""
        return sum(len(self.nested_ideals[k]) for k in range(i)) + j
    
    def __repr__(self):
        subscript_numbers = {0: '₀', 1: '₁', 2: '₂', 3: '₃', 4: '₄', 5: '₅', 6: '₆', 7: '₇', 8: '₈', 9: '₉'}    # TODO: put this in a more global variable
        unk = '�'
        get_subscript = lambda n: subscript_numbers[n] if 0 <= n <= 9 else unk
        
        s_base = self.base.__repr__('I')
        s_ring = f"Ideals of {self.base.R}"
        s_ideals = []
        for i, row in enumerate(self.rows()):
            for j, I in enumerate(row):
                s_ideals.append(f"I{get_subscript(i+1)}{get_subscript(j+1)} = {I}")
        return s_base + "\n" + s_ring + "\n" + "\n".join(s_ideals)


class DoubleNestedHilbertScheme:
    """
    Represents the double nested hilbert scheme corresponding to a given Young diagram.
    The diagram is given in the form of a list of row sizes, 
    for instance the list [3, 3, 2, 1, 1] corresponds to the diagram below.
        ┌─┬─┬─┐
        ├─┼─┼─┤
        ├─┼─┼─┘
        ├─┼─┘
        ├─┤
        └─┘
    """
    
    def __init__(self, diagram: List[int], R: PolyRing = None):
        """
        diagram: a representation of the Young diagram that defines the structure of the space
        """
        # TODO: add sanity check on diagram shape
        if R is None:
            R = infer_poly_ring()
        self.R = R
        self.diagram = diagram
    
    def tangent_space(self, nested_ideals: Union[YoungDiagramIdeals, List[List[PolyRingIdeal]]]):
        if isinstance(nested_ideals, YoungDiagramIdeals):
            return DoubleNestedHilbertSchemeTangentSpace(self, nested_ideals)
        if isinstance(nested_ideals, list):
            return DoubleNestedHilbertSchemeTangentSpace.from_ideal_list(self, nested_ideals)
        raise TypeError
    
    @staticmethod
    def corner_repr(right: bool = False, left: bool = False, top: bool = False, bottom: bool = False):
        if right:
            if top: return '┐'
            if bottom: return '┘'
            return '┤'
        if left:
            if top: return '┌'
            if bottom: return '└'
            return '├'
        if top: return '┬'
        if bottom: return '┴'
        return '┼'

    def row_repr(self, ind: int, length: int, next_length: int = -1, horizontal_pad: int = 5, center_symbol: str = 'Z') -> str:
        subscript_numbers = {0: '₀', 1: '₁', 2: '₂', 3: '₃', 4: '₄', 5: '₅', 6: '₆', 7: '₇', 8: '₈', 9: '₉'}    # TODO: put this in a more global variable
        unk = '�'
        get_subscript = lambda n: subscript_numbers[n] if 0 <= n <= 9 else unk
        
        line = self.corner_repr(left=True, bottom=(0>=next_length))
        middle = '│'
        l = max(length, next_length)
        for i in range(length):
            line += '─'*horizontal_pad
            line += self.corner_repr(right=(i==l-1), bottom=(i>=next_length))
            middle += f" {center_symbol}{get_subscript(ind+1)}{get_subscript(i+1)} │"
        for i in range(length, next_length):
            line += '─'*horizontal_pad
            line += self.corner_repr(right=(i==l-1), top=True)
        return "\n".join((middle, line))

    def __repr__(self, center_symbol: str = 'Z'):
        if len(self.diagram) == 0:
            return ""
        
        next_length = lambda i: (self.diagram[i+1] if i < len(self.diagram)-1 else -1)
        horizontal_pad = 5
        top = self.corner_repr(left=True, top=True)
        top_length = self.diagram[0]
        for i in range(top_length):
            top += '─'*horizontal_pad
            top += self.corner_repr(right=(i==top_length-1), top=True)
        top += "\n"
        s = "\n".join(
            [
                self.row_repr(i, self.diagram[i], next_length(i), horizontal_pad=horizontal_pad, center_symbol=center_symbol)
                for i in range(len(self.diagram))
            ]
        )
        return top + s
    
    def point_from_input(self) -> YoungDiagramIdeals:
        subscript_numbers = {0: '₀', 1: '₁', 2: '₂', 3: '₃', 4: '₄', 5: '₅', 6: '₆', 7: '₇', 8: '₈', 9: '₉'}    # TODO: put this in a more global variable
        unk = '�'
        get_subscript = lambda n: subscript_numbers[n] if 0 <= n <= 9 else unk
        set_global_scope(globals())
        self.R._make_symbols_global_vars(verbose=False)
        
        print(self)
        print(f"Input ideals of {self.R} corresponding to the Zᵢⱼ")

        nested_ideals: List[List[PolyRingIdeal]] = []
        for i, row_length in enumerate(self.diagram):
            nested_ideals.append([])
            for j in range(row_length):
                gens = input(f"Generators of I{get_subscript(i+1)}{get_subscript(j+1)}: ")
                I = eval(f"self.R.ideal({gens})")
                nested_ideals[-1].append(I)
        
        revert_global_scope()
        return YoungDiagramIdeals(self, nested_ideals, self.R)
    
    def point_from_ideal_list(self, nested_ideals: List[List[PolyRingIdeal]]) -> YoungDiagramIdeals:
        diagram_ideals = YoungDiagramIdeals(self, nested_ideals, self.R)
        return diagram_ideals


class DoubleNestedHilbertSchemeTangentSpace:
    def __init__(self, base: DoubleNestedHilbertScheme, diagram_ideals: YoungDiagramIdeals):
        """
        The shape of the nested_ideals list must coincide with the space's diagram
        """
        if base.diagram != diagram_ideals.diagram:
            raise ValueError
        self.base = base
        self.diagram_ideals = diagram_ideals
        self.constraint_sizes: List[int] = []
        self.constraints = self._compute_constraints()
    
    @classmethod
    def from_ideal_list(cls, base: DoubleNestedHilbertScheme, nested_ideals: List[List[PolyRingIdeal]]):
        diagram_ideals = YoungDiagramIdeals(base, nested_ideals, base.R)
        return cls(base, diagram_ideals)
    
    @staticmethod
    def _nested_hom_constraints(I1: PolyRingIdeal = None, I2: PolyRingIdeal = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Required: I₂ subset of I₁
        Constraints for the morphisms I₁ -> O₁ and I₂ -> O₂ to respect the inclusion Z₁ -> Z₂,
        that is commutativity of the following diagram
        I₁ → I₂
        ↓    ↓
        O₁ → O₂
        
        Returns the constraint in two separate matrices, to be glued later
        """
        nested_modules = get_nested_modules(I1, I2)
        S, J1, J2, O1, O2 = nested_modules.components()
        k, m1, m2, n1, n2 = nested_modules.dims()
        phi = HomSpace(J2, J1).get_matrix_representation(lambda f: f)
        psi = HomSpace(O2, O1).get_matrix_representation(lambda f: f)
        C1 = np.zeros((n1*m2, n1*m1))
        C2 = np.zeros((n1*m2, n2*m2))
        for i in range(n1):
            for j in range(m2):
                ind = i*m2 + j
                C1[ind, i*m1: i*m1+m1] += phi[:, j]
                C2[ind, j: n2*m2 + j: m2] -= psi[i]
        return C1, C2

    def _zeros_before(self, ind):
        return sum(self.constraint_sizes[:ind])
    
    def _zeros_after(self, ind):
        return sum(self.constraint_sizes[ind+1:])
    
    def _zeros_between(self, ind1, ind2):
        return sum(self.constraint_sizes[ind1+1:ind2])
    
    def _morphism_constraints(self) -> np.ndarray:
        """
        Constraints for every map I -> O to be a morphism of R-modules.
        Also updates the constraint_sizes attribute
        """
        morphism_constraints = []
        for I in self.diagram_ideals:
            C = HilbertScheme(self.base.R).tangent_space(I).constraints()
            morphism_constraints.append(C)
            size = C.shape[1]
            self.constraint_sizes.append(size)
        
        # padding with zeros
        for ind, C in enumerate(morphism_constraints):
            constraint_count = C.shape[0]
            padding = lambda dim: np.zeros((constraint_count, dim))
            morphism_constraints[ind] = np.concatenate([padding(self._zeros_before(ind)), C, padding(self._zeros_after(ind))], axis=1)
        C = np.concatenate(morphism_constraints, axis=0)
        return C
    
    def _kernel_constraint(self, pos1: Tuple[int, int], pos2: Tuple[int, int]) -> np.ndarray:
        """
        Constraint for the morphisms corresponding I₁ -> O₁ and I₂ -> O₂ to respect the inclusion Z₁ -> Z₂
        Basically a wrapper around _nested_hom_constraints glueing the constraints appropriately
        """
        get_ideal_ind = lambda i, j: (self.diagram_ideals.index_mapping(i, j), self.diagram_ideals[i, j])
        ind1, I1 = get_ideal_ind(*pos1)
        ind2, I2 = get_ideal_ind(*pos2)
        C1, C2 = self._nested_hom_constraints(I1, I2)
        constraint_count = C1.shape[0]
        
        # padding with zeros
        padding = lambda dim: np.zeros((constraint_count, dim))
        C = np.concatenate(
            [
                padding(self._zeros_before(ind1)),
                C1,
                padding(self._zeros_between(ind1, ind2)),
                C2,
                padding(self._zeros_after(ind2))
            ],
            axis=1
        )
        filter_zero = lambda C: C[np.any(C, axis=1)]    # filter out some useless constraints
        C = filter_zero(C)
        return C

    def _compute_constraints(self) -> np.ndarray:
        morphism_constraints = self._morphism_constraints()

        row_constraints = []
        for i, row in enumerate(self.diagram_ideals.rows()):
            for j in range(len(row)-1):
                C = self._kernel_constraint((i, j), (i, j+1))
                row_constraints.append(C)
        
        column_constraints = []
        for j, column in enumerate(self.diagram_ideals.columns()):
            for i in range(len(column)-1):
                C = self._kernel_constraint((i, j), (i+1, j))
                column_constraints.append(C)
        
        constraints = [morphism_constraints] + row_constraints + column_constraints
        C = np.concatenate(
            constraints,
            axis = 0
        )
        return C
    
    def dim(self) -> int:
        max_rank = sum(self.constraint_sizes)
        return max_rank - (np.linalg.matrix_rank(self.constraints) if self.constraints.size > 0 else 0)
    
    def basis(self) -> np.ndarray:
        return null_space(self.constraints)
