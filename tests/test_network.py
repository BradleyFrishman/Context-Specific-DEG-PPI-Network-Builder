"""Check graph semantics against an independent component traversal."""
import random
import unittest
from src.network import networks


def independent_metrics(genes, partners, allowed, combined):
    adj = {g: set() for g in genes}
    for i, a in enumerate(genes):
        for b in genes[i + 1:]:
            shared = bool(partners[a] & partners[b] & allowed)
            direct = b in partners[a] or a in partners[b]
            if shared or (combined and direct):
                adj[a].add(b)
                adj[b].add(a)
    seen, best, edge_count = set(), 0, 0
    for start in genes:
        if start in seen:
            continue
        component, todo = set(), [start]
        while todo:
            g = todo.pop()
            if g in component:
                continue
            component.add(g)
            todo.extend(adj[g] - component)
        seen |= component
        if len(component) > best:
            best = len(component)
            edge_count = sum(len(adj[g] & component) for g in component) // 2
    return best, edge_count


class GraphTests(unittest.TestCase):
    def test_random_graphs_against_independent_traversal(self):
        rng = random.Random(730)
        universe = [f'g{i}' for i in range(20)]
        for _ in range(40):
            partners = {g: {p for p in universe if p != g and rng.random() < .12} for g in universe}
            genes = universe[:10]
            allowed = set(rng.sample(universe, 6))
            masks = {g: sum(1 << universe.index(p) for p in ps & allowed) for g, ps in partners.items()}
            result = networks(genes, partners, masks)
            self.assertEqual(result['second_only'], independent_metrics(genes, partners, allowed, False))
            self.assertEqual(result['first_and_second'], independent_metrics(genes, partners, allowed, True))

    def test_shared_bridge_is_not_an_added_node(self):
        partners = {'a': {'bridge'}, 'b': {'bridge'}, 'c': set()}
        masks = {'a': 1, 'b': 1, 'c': 0}
        self.assertEqual(networks(['a', 'b', 'c'], partners, masks)['second_only'], (2, 1))

    def test_empty_graph(self):
        self.assertEqual(networks([], {}, {})['second_only'], (0, 0))


if __name__ == '__main__':
    unittest.main()
