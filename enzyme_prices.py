"""Public manufacturer prices for exact restriction-enzyme products."""
import gzip
import json
import re
from dataclasses import dataclass
from decimal import Decimal
from urllib.parse import urlencode, urlsplit
from urllib.request import Request, urlopen

from genome_download import HttpTransport

SIBENZYME = 'SibEnzyme'
CLONING = 'Клонинг Фасилити'
SUPPLIERS = (SIBENZYME, CLONING)
_SIB_API = 'https://sibenzyme.com/wp-json/wc/store/v1/products'
_CLONING_BASE = 'https://cloning.tech/products/reagents/'
_CLONING_PATHS = {
    'BamHI': 'bamhi-restrictase', 'DpnI': 'dpni-restrictase', 'DraI': 'drai-restrictase',
    'EcoRI': 'ecori-restrictase', 'HgaI': 'hgai-restrictase', 'HindIII': 'hindiii-restrictase',
    'NotI': 'noti-cutfast', 'XbaI': 'xbai-restrictase',
}
_JSON_LD = re.compile(r'<script\b[^>]*\btype=["\']application/ld\+json["\'][^>]*>(.*?)</script>',
                      re.IGNORECASE | re.DOTALL)


@dataclass(frozen=True)
class Quote:
    sku: str
    pack: str
    amount: Decimal
    currency: str
    url: str


def sources_for(enzyme):
    return SUPPLIERS if enzyme in _CLONING_PATHS else (SIBENZYME,)


def is_product_link(supplier, url):
    parts = urlsplit(url)
    host, prefix = (('sibenzyme.com', '/product/') if supplier == SIBENZYME
                    else ('cloning.tech', '/products/reagents/'))
    return parts.scheme == 'https' and parts.netloc == host and parts.path.startswith(prefix)


def _read(url):
    request = Request(url, headers={'User-Agent': 'RFLP-Picker', 'Accept-Encoding': 'gzip',
                                    'Cache-Control': 'no-cache'})
    with urlopen(request, timeout=8, context=HttpTransport._ssl_context()) as response:
        stream = gzip.GzipFile(fileobj=response) if response.headers.get('Content-Encoding') == 'gzip' else response
        data = stream.read(4_000_001)
    if len(data) > 4_000_000:
        raise ValueError('Manufacturer response is too large')
    return data.decode('utf-8')


def _sibenzyme(enzyme, read):
    # The catalog inserts spaces before the Roman numeral, e.g. "EcoR I".
    stem = re.sub(r'[IVX]+$', '', enzyme) or enzyme
    products = json.loads(read(f'{_SIB_API}?{urlencode({"search": stem, "per_page": 50})}'))
    matches = [p for p in products if ''.join(p['name'].split()).casefold() == enzyme.casefold()]
    if not matches:
        return ()
    if len(matches) != 1:
        raise ValueError('Ambiguous SibEnzyme product')
    product, = matches
    url = product['permalink']
    if not is_product_link(SIBENZYME, url):
        raise ValueError('Unexpected SibEnzyme product link')
    if product['type'] == 'variable':
        params = urlencode({'parent': product['id'], 'type': 'variation', 'per_page': 100})
        variants = json.loads(read(f'{_SIB_API}?{params}'))
        if not variants or any(v['parent'] != product['id'] for v in variants):
            raise ValueError('Unexpected SibEnzyme variations')
    else:
        variants = [product]
    quotes = []
    for variant in variants:
        price = variant['prices']
        if not variant['sku'] or not price['price']:
            continue
        attributes = dict(pair.split(': ', 1) for pair in variant['variation'].split(', ') if ': ' in pair)
        pack = ' · '.join(attributes[key] for key in ('Фасовка', 'Вариант') if key in attributes)
        amount = Decimal(price['price']).scaleb(-price['currency_minor_unit'])
        quotes.append(Quote(variant['sku'], pack, amount, price['currency_code'], url))
    return tuple(sorted(quotes, key=lambda quote: (quote.amount, quote.sku)))


def _cloning(enzyme, read):
    url = f'{_CLONING_BASE}{_CLONING_PATHS[enzyme]}/'
    products = []
    for script in _JSON_LD.findall(read(url)):
        data = json.loads(script)
        nodes = data.get('@graph', [data]) if isinstance(data, dict) else data
        products.extend(node for node in nodes if isinstance(node, dict) and node.get('@type') == 'Product')
    matches = [product for product in products if product['name'].split()[1:2] == [enzyme]]
    if len(matches) != 1:
        raise ValueError('Expected exactly one Cloning Facility product for this enzyme')
    product, = matches
    quotes = []
    for offer in product['offers']:
        pack, _, catalog_number = offer['name'].partition(', Кат.№ ')
        quotes.append(Quote(catalog_number or offer['sku'], pack, Decimal(offer['price']),
                            offer['priceCurrency'], url))
    return tuple(quotes)


def fetch_prices(supplier, enzyme, read=_read):
    """Populate only exact products. Transport and parse errors propagate to the caller."""
    fetch = {SIBENZYME: _sibenzyme, CLONING: _cloning}[supplier]
    quotes = fetch(enzyme, read)
    if any(not quote.amount.is_finite() or quote.amount <= 0 or not quote.currency for quote in quotes):
        raise ValueError('Invalid manufacturer price')
    return quotes
