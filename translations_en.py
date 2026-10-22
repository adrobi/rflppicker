"""English texts for the interface (i18n spec §3–4).

EN maps each Russian source passed to t() or N_() to its English text, with
the same {placeholders}. CORE lists (regular expression over the core's Russian
text, English replacement) pairs; core_text() applies them in order, so longer
phrases come before the shorter ones they contain.
"""

EN = {
}

CORE = [
    # Reasons in the result rows.
    (r'SNP ближе (\d+) п\.н\. к краю окна; увеличьте фланг',
     r'SNP is closer than \1 bp to the window edge; increase the flank'),
    (r'Primer3 не вернул подходящих пар', 'Primer3 returned no suitable pairs'),
    (r'Не найдено сочетание для различения трёх генотипов', 'No primer pair and enzyme tell the three genotypes apart'),
    (r'Подходящий фермент не найден', 'No suitable enzyme found'),
    (r'; окон с неопределёнными основаниями: (\d+)', r'; windows with ambiguous bases: \1'),
    (r'; ошибок расчёта ферментов: (\d+)', r'; enzyme calculation errors: \1'),
    (r'Три генотипа различимы во всём заданном диапазоне разрешения',
     'The three genotypes are distinguishable over the whole resolution range'),
    (r'Три генотипа различимы только при исходном разрешении',
     'The three genotypes are distinguishable at the nominal resolution only'),
    (r'Три генотипа неразличимы по присутствию полос', 'The three genotypes cannot be told apart by band presence'),
    (r'Оценка окна без подобранных праймеров\. ', 'Window estimate without designed primers. '),
    (r'dCAPS: несовпадение в праймере создаёт сайт в аллеле (REF|ALT)\. ',
     r'dCAPS: the primer mismatch creates the site in the \1 allele. '),
    (r'Модуль Primer3 недоступен', 'The Primer3 module is not available'),
    (r'Primer3 вернул праймеры вне окна или перекрывающие SNP', 'Primer3 returned primers outside the window or over the SNP'),
    (r'Несогласованные размеры продукта/праймеров Primer3', 'Inconsistent Primer3 product or primer sizes'),
    (r'Последовательности праймеров не совпадают с референсом', 'Primer sequences do not match the reference'),
    (r'Primer3 вернул некорректную оценку', 'Primer3 returned an invalid score'),
    (r'REF не совпадает с FASTA', 'REF does not match the FASTA'),
    # Details of rejected lines.
    (r'Ожидалось минимум 4 поля', 'At least 4 fields expected'),
    (r'Поддерживаются только одиночные SNP A/C/G/T', 'Only single-base A/C/G/T SNPs are supported'),
    (r"CHROM '([^']*)' отсутствует в FASTA", r"CHROM '\1' is not in the FASTA"),
    (r'Позиция вне FASTA', 'Position outside the FASTA'),
    # Parameter checks (validate_params) and the band model.
    (r'Не задан параметр: ', 'Parameter not set: '),
    (r'Файл не найден \((\w+)\): ', r'File not found (\1): '),
    (r'Некорректный параметр: ', 'Invalid parameter: '),
    (r'Минимальный размер продукта больше максимального', 'Minimum product size is greater than maximum'),
    (r'Минимальный размер фрагмента больше максимального', 'Minimum fragment size is greater than maximum'),
    (r'Минимальный Tm больше максимального', 'Minimum Tm is greater than maximum'),
    (r': требуется целое число от (-?\d+) до (-?\d+)', r': an integer from \1 to \2 is required'),
    (r'dcaps: допустимо off, fallback или always', 'dcaps: allowed values are off, fallback or always'),
    (r'Параметры модели полос должны быть конечными числами', 'Band model parameters must be finite numbers'),
    (r'Некорректный диапазон видимых фрагментов', 'Invalid range of visible fragments'),
    (r'Некорректное разрешение геля', 'Invalid gel resolution'),
    (r'Коэффициент ухудшения разрешения должен быть от 1 до 5', 'The resolution stress factor must be between 1 and 5'),
    (r'Длины фрагментов должны быть положительными', 'Fragment lengths must be positive'),
    # The run log (status_cb).
    (r'Поиск: (\S+); SNP: (\d+); ферментов: (\d+)', r'Search: \1; SNPs: \2; enzymes: \3'),
    (r'Оценка полос: модель присутствия, полная рестрикция; '
     r'яркость и специфичность ПЦР по всему геному не оцениваются\.',
     'Band model: presence of bands, complete digestion; '
     'band brightness and genome-wide PCR specificity are not assessed.'),
    (r'\[INFO\] Отбор: ', '[INFO] Filters: '),
    (r'только продаваемые ферменты', 'commercial enzymes only'),
    (r'изошизомеры объединены \(одна строка на сайт\)', 'isoschizomers grouped (one row per site)'),
    (r'без сайта фермента в праймерах', 'no enzyme site in primers'),
    (r'повторы \(строчные буквы FASTA\) не на 3′-концах праймеров', 'repeats (lower-case FASTA) not at primer 3′ ends'),
    (r'без фоновых вариантов, меняющих разрезы', 'no background variants that change cuts'),
    (r'SNP не ближе (\d+) п\.н\. к концам ампликона', r'SNP at least \1 bp from amplicon ends'),
    (r'dCAPS если нет robust-схемы, несовпадений до (\d+)', r'dCAPS if no robust assay, mismatches up to \1'),
    (r'dCAPS для каждого SNP, несовпадений до (\d+)', r'dCAPS for every SNP, mismatches up to \1'),
    (r'Фоновые варианты: (\d+) в окнах SNP из (\d+) записей VCF',
     r'Background variants: \1 in SNP windows out of \2 VCF records'),
    (r': пар (\d+), кандидатов (\d+), неразличимых (\d+)', r': pairs \1, candidates \2, indistinguishable \3'),
    (r'В окнах SNP нет строчных букв: FASTA без soft-mask, маскирование повторов не действует\.',
     'No lower-case letters in SNP windows: the FASTA is not soft-masked, repeat masking has no effect.'),
    (r'NORMALIZE: загружено вариантов: ', 'NORMALIZE: variants loaded: '),
]
