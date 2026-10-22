"""English texts for the interface (i18n spec §3–4).

EN maps each Russian source passed to t() or N_() to its English text, with
the same {placeholders}. CORE lists (regular expression over the core's Russian
text, English replacement) pairs; core_text() applies them in order, so longer
phrases come before the shorter ones they contain.
"""

EN = {
    # snp_summary: statuses, notes, texts for the lab notebook, Excel sheets.
    'Надёжно': 'Robust',
    'Только на хорошем геле': 'Good gel only',
    'Нет схемы': 'No assay',
    'Ошибка данных': 'Data error',
    'Хромосома не найдена в геноме': 'Chromosome not found in the genome',
    'Позиция за пределами хромосомы': 'Position beyond the chromosome end',
    'Не одиночная замена A/C/G/T': 'Not a single A/C/G/T substitution',
    'REF похож на обратную цепь': 'REF looks like the reverse strand',
    'ALT совпадает с геномом': 'ALT matches the genome',
    'REF не совпадает с геномом': 'REF does not match the genome',
    'Попробуйте расширить диапазон продукта или Tm': 'Try a wider product size range or Tm',
    'Проверьте параметры Primer3 в окне «Дополнительно»': 'Check the Primer3 parameters in «Advanced»',
    'Строка не прочитана: {details}': 'Line not read: {details}',
    'REF и ALT переставлены по геному: REF = {ref}': 'REF and ALT swapped to match the genome: REF = {ref}',
    'Три генотипа неразличимы по присутствию полос': 'The three genotypes cannot be told apart by band presence',
    'Выбрано вручную, лучшая — {enzyme}': 'Chosen by hand, best is {enzyme}',
    'Только при исходном разрешении; труднее всего: {pair}': 'Nominal resolution only; hardest pair: {pair}',
    'Сайт фермента в праймере': 'Enzyme site in a primer',
    'Фоновые варианты меняют разрезы: {variants}': 'Background variants change the cuts: {variants}',
    'Известных вариантов в ампликоне: {count}': 'Known variants in the amplicon: {count}',
    'Повторы {percent}% ампликона': 'Repeats cover {percent}% of the amplicon',
    'Поиск по окну, праймеры не подобраны': 'Window search, no primers designed',
    'Нет результата для этой строки': 'No result for this line',
    'Попробуйте другой диапазон продукта': 'Try another product size range',
    'Попробуйте dCAPS «Для каждого SNP»': 'Try dCAPS «For every SNP»',
    'Статус: {status}': 'Status: {status}',
    'Фермент: {enzyme}, сайт {site}': 'Enzyme: {enzyme}, site {site}',
    'Продукт: {size} п.н.': 'Product: {size} bp',
    'Полосы: REF/REF {ref_ref}; REF/ALT {ref_alt}; ALT/ALT {alt_alt}':
        'Bands: REF/REF {ref_ref}; REF/ALT {ref_alt}; ALT/ALT {alt_alt}',
    'Примечание: {note}': 'Note: {note}',
    'Труднее всего различить: {pair}': 'Hardest to tell apart: {pair}',
    '; оценка при ухудшении {score}': '; stressed score {score}',
    'Хромосома в геноме': 'Genome chromosome',
    'Статус': 'Status',
    'Выбор': 'Choice',
    'Фермент': 'Enzyme',
    'Сайт': 'Site',
    'Изошизомеры': 'Isoschizomers',
    'Поставщики': 'Suppliers',
    'Тип схемы': 'Assay type',
    'Праймер F 5′→3′': 'Primer F 5′→3′',
    'Праймер R 5′→3′': 'Primer R 5′→3′',
    'Продукт, п.н.': 'Product, bp',
    'Полосы REF/REF': 'Bands REF/REF',
    'Полосы REF/ALT': 'Bands REF/ALT',
    'Полосы ALT/ALT': 'Bands ALT/ALT',
    'Оценка': 'Score',
    'Оценка при ухудшении': 'Stressed score',
    'Примечание': 'Note',
    'Строка файла': 'File line',
    'Имя': 'Name',
    'Последовательность 5′→3′': 'Sequence 5′→3′',
    'Длина': 'Length',
    'природный сайт': 'natural site',
    'вручную': 'by hand',
    'авто': 'auto',
    'Сводка': 'Summary',
    'Праймеры': 'Primers',
    # result_model: titles of the core columns.
    'Ранг для SNP': 'Rank for the SNP',
    'Пара Primer3 №': 'Primer3 pair no.',
    'Штраф Primer3': 'Primer3 penalty',
    'Проверено пар': 'Pairs checked',
    'Различимость генотипов': 'Genotype separation',
    'Оценка различимости (> 1)': 'Separation score (> 1)',
    'Оценка при ухудшении (> 1)': 'Stressed score (> 1)',
    'Труднее всего различить': 'Hardest to tell apart',
    'Полосы REF/REF при ухудшении': 'Bands REF/REF, stressed',
    'Полосы REF/ALT при ухудшении': 'Bands REF/ALT, stressed',
    'Полосы ALT/ALT при ухудшении': 'Bands ALT/ALT, stressed',
    'Число полос в гетерозиготе': 'Bands in the heterozygote',
    'Невидимых фрагментов REF': 'Invisible REF fragments',
    'Невидимых фрагментов ALT': 'Invisible ALT fragments',
    'Режим поиска': 'Search mode',
    'Модель оценки': 'Scoring model',
    'Вариант': 'Variant',
    'Контиг': 'Contig',
    'Паттерн': 'Pattern',
    'Фрагменты REF, п.н.': 'REF fragments, bp',
    'Фрагменты ALT, п.н.': 'ALT fragments, bp',
    'Различие, п.н.': 'Difference, bp',
    'Начало (1-based)': 'Start (1-based)',
    'Конец (1-based)': 'End (1-based)',
    'Длина, п.н.': 'Length, bp',
    'Левый праймер 5′→3′': 'Left primer 5′→3′',
    'Правый праймер 5′→3′': 'Right primer 5′→3′',
    'Tm лев., °C': 'Tm left, °C',
    'Tm прав., °C': 'Tm right, °C',
    'Продукт Primer3': 'Primer3 product',
    'Начало L': 'Start L',
    'Длина L': 'Length L',
    'Начало R': 'Start R',
    'Длина R': 'Length R',
    'Пояснение': 'Explanation',
    'Порог Δ, п.н.': 'Δ threshold, bp',
    'Сайт в праймере': 'Site in a primer',
    'Повторы, %': 'Repeats, %',
    'Известных вариантов в ампликоне': 'Known variants in the amplicon',
    'Варианты, меняющие разрезы/полосы': 'Variants that change cuts or bands',
    'Несовпадения в праймере': 'Primer mismatches',
    'Строка': 'Line',
    'Хромосома': 'Chromosome',
    'Позиция': 'Position',
    'Подробности': 'Details',
    'Исходная строка': 'Original line',
    # rflp_cli: help and messages.
    '''RFLP Picker из командной строки: тот же анализ, что в окне программы.

Пример:
    python rflp_cli.py --fasta genome.fa --variants snps.tsv --out results.xlsx

Результат сохраняется в CSV или Excel (по расширению --out). В Excel, как и
при экспорте из окна, добавляются листы Run (параметры, версии, сводка) и
Rejected (отклонённые варианты). Журнал выводится в stderr.''':
    '''RFLP Picker from the command line: the same analysis as in the program window.

Example:
    python rflp_cli.py --fasta genome.fa --variants snps.tsv --out results.xlsx

The result is saved as CSV or Excel (by the extension of --out). Like the
export from the window, Excel gets the Run sheet (parameters, versions, totals)
and the Rejected sheet (rejected variants). The log goes to stderr.''',
    'Язык справки, сообщений и листов Сводка/Праймеры (по умолчанию — язык системы)':
        'Language of the help, the messages and the Summary/Primers sheets (default: the system language)',
    'данные': 'data',
    'Референсный геном FASTA (рядом создаётся .fai)': 'Reference genome FASTA (a .fai index is created next to it)',
    'CSV/TSV/TXT с колонками chrom, pos, ref, alt': 'CSV/TSV/TXT with the columns chrom, pos, ref, alt',
    'Файл соответствий названий хромосом': 'File that maps chromosome names',
    'Файл результата: .csv или .xlsx': 'Result file: .csv or .xlsx',
    'CSV с отклонёнными вариантами (в .xlsx они всегда на листе Rejected)':
        'CSV of rejected variants (in .xlsx they are always on the Rejected sheet)',
    'Исходный набор параметров: recommended — настройки по итогам '
    'проверки (как в окне программы), v1.1 — воспроизведение v1.1.0 '
    '(по умолчанию: %(default)s)':
        'Starting parameter set: recommended — the validated settings (as in the program window), '
        'v1.1 — reproduces v1.1.0 (default: %(default)s)',
    'поиск (по умолчанию — значения пресета)': 'search (default: the preset values)',
    'Фланг вокруг SNP, п.н.': 'Flank around the SNP, bp',
    'Мин. видимый фрагмент, п.н.': 'Min. visible fragment, bp',
    'Макс. видимый фрагмент, п.н.': 'Max. visible fragment, bp',
    'Диагностическая дельта, п.н.': 'Diagnostic delta, bp',
    'Макс. число разрезов в ампликоне': 'Max. number of cuts in the amplicon',
    'Только появление/исчезновение сайта': 'Site gain or loss only',
    'Поставщики через запятую: коды REBASE (N,E) или названия':
        'Suppliers, comma-separated: REBASE codes (N,E) or names',
    'отбор и проверки': 'filters and checks',
    'Только ферменты, у которых в REBASE есть поставщик': 'Only enzymes with a supplier in REBASE',
    'Одна строка на сайт узнавания; остальные ферменты — в столбце isoschizomers':
        'One row per recognition site; the other enzymes go to the isoschizomers column',
    'SNP не ближе BP п.н. к каждому концу ампликона (0 — не ограничивать)':
        'SNP at least BP bp from each amplicon end (0: no limit)',
    'Исключать схемы, где праймер сам содержит сайт фермента':
        'Exclude assays where a primer itself contains the enzyme site',
    'Повторы (строчные буквы soft-mask FASTA) не на 3′-концах праймеров':
        'Repeats (lower-case soft-masked FASTA) not at primer 3′ ends',
    'VCF/VCF.gz известных вариантов: праймеры их избегают, влияние на разрезы '
    'отмечается (названия хромосом — через --names)':
        'VCF/VCF.gz of known variants: primers avoid them, effects on cuts are '
        'flagged (chromosome names via --names)',
    'Исключать схемы, где фоновые варианты меняют разрезы или длины полос':
        'Exclude assays where background variants change cuts or band sizes',
    'Праймеры с несовпадением (dCAPS): fallback — только если нет robust-схемы':
        'Mismatch primers (dCAPS): fallback — only when there is no robust assay',
    'Наибольшее число несовпадений в праймере dCAPS': 'Maximum number of mismatches in a dCAPS primer',
    'праймеры и совместный подбор': 'primers and joint design',
    'Подбирать праймеры Primer3 (--no-primer3 — поиск по окну)':
        'Design primers with Primer3 (--no-primer3: window search)',
    'Совместный подбор праймеров и фермента': 'Choose primers and enzyme together',
    'Пар Primer3 на SNP': 'Primer3 pairs per SNP',
    'Результатов на SNP (0 — все)': 'Results per SNP (0: all)',
    'Мин. размер продукта, п.н.': 'Min. product size, bp',
    'Макс. размер продукта, п.н.': 'Max. product size, bp',
    'Мин. Tm, °C': 'Min. Tm, °C',
    'Макс. Tm, °C': 'Max. Tm, °C',
    'модель геля': 'gel model',
    'Разрешение, п.н.': 'Resolution, bp',
    'Разрешение, %%': 'Resolution, %%',
    'Ухудшение разрешения, ×': 'Resolution stress, ×',
    'нормализация': 'normalization',
    'Не исправлять REF/ALT (проверка REF по FASTA остаётся)': 'Do not fix REF/ALT (REF is still checked against the FASTA)',
    'Сохранять cleaned/rejected в results/normalize': 'Save cleaned/rejected to results/normalize',
    'Не выводить журнал': 'Do not print the log',
    '--out: требуется расширение .csv или .xlsx': '--out: the extension must be .csv or .xlsx',
    'Ошибка параметров: {error}': 'Parameter error: {error}',
    'Анализ прерван.': 'Analysis interrupted.',
    '[INFO] Сохранено: {out} · строк {rows} · SNP принято {accepted}, с кандидатами {candidates}, '
    'с robust-кандидатами {robust}':
        '[INFO] Saved: {out} · rows {rows} · SNPs accepted {accepted}, with candidates {candidates}, '
        'with robust candidates {robust}',
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
