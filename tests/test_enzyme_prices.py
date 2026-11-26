import json
import unittest
from decimal import Decimal
from urllib.parse import parse_qs, urlsplit

from enzyme_prices import CLONING, SIBENZYME, fetch_prices


class PriceTests(unittest.TestCase):
    def test_sibenzyme_matches_exact_enzyme_and_prices_each_variant(self):
        card = 'https://sibenzyme.com/product/acs-i/'
        payloads = {
            'search': [
                {'name': 'Acs I', 'id': 42, 'permalink': card, 'type': 'variable'},
                {'name': 'SE-буфер AcsI', 'id': 43, 'type': 'simple'},
                {'name': 'Aci I', 'id': 44, 'type': 'variable'},
            ],
            'parent': [
                {'parent': 42, 'sku': 'E013', 'variation': 'Фасовка: 500 е.а., Вариант: Regular',
                 'prices': {'price': '650000', 'currency_minor_unit': 2, 'currency_code': 'RUB'}},
                {'parent': 42, 'sku': 'E014', 'variation': 'Фасовка: 2500 е.а.',
                 'prices': {'price': '2600000', 'currency_minor_unit': 2, 'currency_code': 'RUB'}},
            ],
        }

        def read(url):
            query = parse_qs(urlsplit(url).query)
            if 'parent' in query:
                self.assertEqual(query['parent'], ['42'])
                return json.dumps(payloads['parent'])
            return json.dumps(payloads['search'])

        quotes = fetch_prices(SIBENZYME, 'AcsI', read)
        self.assertEqual([(q.sku, q.pack, q.amount, q.currency, q.url) for q in quotes], [
            ('E013', '500 е.а. · Regular', Decimal('6500'), 'RUB', card),
            ('E014', '2500 е.а.', Decimal('26000'), 'RUB', card),
        ])

    def test_cloning_facility_reads_product_offers_not_other_json(self):
        payload = {'@type': 'Product', 'name': 'Рестриктаза DraI CutFast', 'offers': [
            {'@type': 'Offer', 'sku': 'drai-50', 'name': '50 реакций, 50 мкл, Кат.№ DRA1S',
             'price': '5000', 'priceCurrency': 'RUB'},
        ]}
        graph = {'@context': 'https://schema.org', '@graph': [{'@type': 'BreadcrumbList'}, payload]}
        html = f'<script type="application/ld+json">{json.dumps(graph)}</script>'
        quote, = fetch_prices(CLONING, 'DraI', lambda url: html)
        self.assertEqual((quote.sku, quote.pack, quote.amount), ('DRA1S', '50 реакций, 50 мкл', Decimal(5000)))
        self.assertEqual(quote.url, 'https://cloning.tech/products/reagents/drai-restrictase/')

    def test_cloning_wrong_product_is_a_source_error(self):
        html = '<script type="application/ld+json">{"@type":"Product","name":"Рестриктаза EcoRI"}</script>'
        with self.assertRaises(ValueError):
            fetch_prices(CLONING, 'DraI', lambda url: html)

    def test_no_exact_sibenzyme_result_is_not_a_quote(self):
        self.assertEqual(fetch_prices(SIBENZYME, 'AcsI', lambda url: '[{"name":"Asc I","id":9}]'), ())


if __name__ == '__main__':
    unittest.main()
