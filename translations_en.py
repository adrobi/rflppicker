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
    'Проверьте параметры Primer3 в окне «Дополнительно»': 'Check the Primer3 parameters in “Advanced”',
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
    'Попробуйте другой диапазон продукта или Tm': 'Try another product range or Tm',
    'Попробуйте dCAPS «Для каждого SNP»': 'Try dCAPS “For every SNP”',
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
Rejected (отклонённые варианты). Журнал выводится в stderr.

Скачать референсный геном с NCBI в папку genoms рядом с программой:
    python rflp_cli.py --download-genome GCF_016772045.1''':
    '''RFLP Picker from the command line: the same analysis as in the program window.

Example:
    python rflp_cli.py --fasta genome.fa --variants snps.tsv --out results.xlsx

The result is saved as CSV or Excel (by the extension of --out). Like the
export from the window, Excel gets the Run sheet (parameters, versions, totals)
and the Rejected sheet (rejected variants). The log goes to stderr.

Download a reference genome from NCBI into the genoms folder next to the program:
    python rflp_cli.py --download-genome GCF_016772045.1''',
    'Язык справки, сообщений и листов Сводка/Праймеры (по умолчанию — язык системы)':
        'Language of the help, the messages and the Summary/Primers sheets (default: the system language)',
    'данные': 'data',
    'Референсный геном FASTA (рядом создаётся .fai)': 'Reference genome FASTA (a .fai index is created next to it)',
    'CSV/TSV/TXT с колонками chrom, pos, ref, alt': 'CSV/TSV/TXT with the columns chrom, pos, ref, alt',
    'Файл соответствий названий хромосом': 'File that maps chromosome names',
    'Файл результата: .csv или .xlsx': 'Result file: .csv or .xlsx',
    'CSV с отклонёнными вариантами (в .xlsx они всегда на листе Rejected)':
        'CSV of rejected variants (in .xlsx they are always on the Rejected sheet)',
    'поиск (по умолчанию — как в окне программы)': 'search (default: as in the program window)',
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
    'dCAPS: верхняя граница Tm обоих праймеров, °C (нижняя — --tm-min)':
        'dCAPS: upper Tm limit of both primers, °C (the lower one is --tm-min)',
    'dCAPS: мин. размер продукта, п.н.': 'dCAPS: min. product size, bp',
    'dCAPS: макс. размер продукта, п.н.': 'dCAPS: max. product size, bp',
    'dCAPS: порог Δ, п.н. (вместо --delta для схем dCAPS)':
        'dCAPS: Δ threshold, bp (instead of --delta for dCAPS assays)',
    'праймеры': 'primers',
    'Подбирать праймеры Primer3 (--no-primer3 — поиск по окну)':
        'Design primers with Primer3 (--no-primer3: window search)',
    'Пар Primer3 на SNP': 'Primer3 pairs per SNP',
    'Результатов на SNP (0 — все)': 'Results per SNP (0: all)',
    'Мин. размер продукта, п.н.': 'Min. product size, bp',
    'Макс. размер продукта, п.н.': 'Max. product size, bp',
    'Мин. Tm, °C': 'Min. Tm, °C',
    'Макс. Tm, °C': 'Max. Tm, °C',
    'модель геля': 'gel model',
    'Разрешение, п.н.': 'Resolution, bp',
    'Разрешение, %%': 'Resolution, %%',
    'Ухудшение разрешения, ×': 'Resolution stress factor, ×',
    'нормализация': 'normalization',
    'Не исправлять REF/ALT (проверка REF по FASTA остаётся)':
        'Do not fix REF/ALT (REF is still checked against the FASTA)',
    'Сохранять cleaned/rejected в results/normalize': 'Save cleaned/rejected to results/normalize',
    'Не выводить журнал': 'Do not print the log',
    '--out: требуется расширение .csv или .xlsx': '--out: the extension must be .csv or .xlsx',
    'Ошибка параметров: {error}': 'Parameter error: {error}',
    'Анализ прерван.': 'Analysis interrupted.',
    '{key}: папка не существует: {folder}': '{key}: the folder does not exist: {folder}',
    'Ошибка: {error}': 'Error: {error}',
    # specificity
    'Исключать пары': 'Exclude pairs', 'Только отмечать': 'Flag only', 'Не проверять': 'Do not check',
    'Специфичность': 'Specificity',
    'Праймеры ищутся по всему геному: пара, которая даёт продукт ПЦР ещё в '
    'одном месте, исключается или получает примечание. Проверка читает '
    'весь геном и добавляет к расчёту около минуты на каждые 3 млрд оснований.':
        'The primers are looked up in the whole genome: a pair that gives a PCR product in one more place '
        'is excluded or gets a note. The check reads the whole genome and adds about a minute per '
        '3 billion bases to the analysis.',
    'Специфичность: несовпадений, до': 'Specificity: mismatches, up to',
    'Место отжига праймера в геноме: 12 оснований 3′-конца совпадают точно, в остальной '
    'части праймера — не больше этого числа несовпадений.':
        'An annealing site of a primer in the genome: the 12 bases of the 3′ end match exactly, the rest '
        'of the primer has at most this number of mismatches.',
    'Специфичность: продукт до, п.н.': 'Specificity: product up to, bp',
    'Посторонний продукт ПЦР: два места отжига, обращённые друг к другу, не дальше '
    'этого расстояния.':
        'Another PCR product: two annealing sites facing each other, not further apart than this.',
    'Поиск праймеров по всему геному: exclude — исключать пары с посторонним продуктом ПЦР, '
    'warn — только отмечать, off — не проверять':
        'Look the primers up in the whole genome: exclude — leave out pairs with another PCR product, '
        'warn — only flag them, off — do not check',
    'Специфичность: несовпадений вне 12 оснований 3′-конца, не больше':
        'Specificity: mismatches outside the 12 bases of the 3′ end, at most',
    'Специфичность: наибольший посторонний продукт, п.н.': 'Specificity: the largest other product, bp',
    'Посторонние продукты ПЦР: {sizes} п.н.': 'Other PCR products: {sizes} bp',
    'Мест отжига левого праймера': 'Annealing sites of the left primer',
    'Мест отжига правого праймера': 'Annealing sites of the right primer',
    'Посторонних продуктов ПЦР': 'Other PCR products',
    'Длины посторонних продуктов, п.н.': 'Sizes of the other products, bp',
    # laboratory sheet
    'Лист PDF': 'PDF sheet',
    'Лабораторный лист выбранной схемы для печати: праймеры, фермент, полосы, гель и карта':
        'A laboratory sheet of the chosen assay for printing: primers, enzyme, bands, gel and map',
    'Листы PDF: все SNP': 'PDF sheets: all SNPs',
    'Листы PDF: только показанные': 'PDF sheets: shown only',
    'Нет SNP со схемой: листы не созданы.': 'No SNP has an assay: no sheets were made.',
    'Дата': 'Date', 'Оператор': 'Operator', '№ геля': 'Gel no.',
    'RFLP Picker · лабораторный лист': 'RFLP Picker · laboratory sheet',
    'Фото геля': 'Gel photo',
    'геном {genome}': 'genome {genome}',
    'продукт {prod_min}–{prod_max} п.н., Tm {tm_min}–{tm_max} °C': 'product {prod_min}–{prod_max} bp, Tm {tm_min}–{tm_max} °C',
    'лист создан {date}': 'sheet made on {date}',
    'Не удалось создать PDF: {path}': 'The PDF could not be made: {path}',
    # genome download
    'Овца': 'Sheep', 'Корова': 'Cattle', 'Коза': 'Goat', 'Свинья': 'Pig', 'Курица': 'Chicken', 'Рис': 'Rice',
    'Пшеница': 'Wheat', 'Ячмень': 'Barley', 'Сорго': 'Sorghum',
    'Скачать геном…': 'Download a genome…',
    'Скачать геном': 'Download a genome',
    'Сборка генома на NCBI': 'Genome assembly at NCBI',
    '{species} — {name} ({accession}, {size} ГБ)': '{species} — {name} ({accession}, {size} GB)',
    'Другая сборка…': 'Another assembly…',
    'Номер сборки, например GCF_016772045.1': 'Assembly accession, for example GCF_016772045.1',
    'Координаты в файле SNP должны относиться к той же сборке. Геном сохраняется в папку '
    'genoms рядом с программой; файл соответствий названий хромосом создаётся автоматически.':
        'The coordinates of the SNP file must belong to the same assembly. The genome is saved to the '
        'genoms folder next to the program; the chromosome name mapping file is made automatically.',
    'Отмена': 'Cancel',
    'Скачать': 'Download',
    'Останавливаю загрузку…': 'Stopping the download…',
    'Загрузка с NCBI: {done} из {total} МБ': 'Downloading from NCBI: {done} of {total} MB',
    'Загрузка остановлена; при повторе она продолжится с этого места.':
        'The download was stopped; a repeat continues from this point.',
    'Загрузка с NCBI…': 'Downloading from NCBI…',
    'Проверка контрольной суммы…': 'Checking the checksum…',
    'Распаковка…': 'Unpacking…',
    'Индексация…': 'Indexing…',
    'NCBI ответил ошибкой {code} ({url})': 'NCBI answered with error {code} ({url})',
    'Нет связи с NCBI: {error}': 'No connection to NCBI: {error}',
    'NCBI не поддержал докачку; удалите незавершённый файл и повторите':
        'NCBI did not continue the download; delete the unfinished file and repeat',
    'Соединение с NCBI прервано: {error}': 'The connection to NCBI was lost: {error}',
    'Номер сборки должен иметь вид GCF_016772045.1, а не «{value}»':
        'The assembly accession must look like GCF_016772045.1, not “{value}”',
    'Сборка {accession} не найдена на NCBI': 'The assembly {accession} was not found at NCBI',
    'Неожиданное имя сборки на NCBI: {name}': 'An unexpected assembly name at NCBI: {name}',
    'Соединение с NCBI прервано: получено {done} из {total} МБ. '
    'Повторите загрузку — она продолжится с этого места.':
        'The connection to NCBI was lost: {done} of {total} MB received. '
        'Repeat the download: it continues from this point.',
    'На NCBI нет контрольной суммы файла {name}': 'NCBI has no checksum of the file {name}',
    'Недостаточно места на диске: нужно около {size} ГБ в папке {folder}':
        'Not enough disk space: about {size} GB are needed in the folder {folder}',
    'Контрольная сумма не совпала с NCBI: загрузка повреждена, повторите её':
        'The checksum differs from NCBI: the download is damaged, repeat it',
    'Геном: {fasta}\nХромосомы: {names}': 'Genome: {fasta}\nChromosomes: {names}',
    '{key}: это папка, а не файл: {path}': '{key}: this is a folder, not a file: {path}',
    '{key}: нет права записи в папку: {folder}': '{key}: no permission to write to the folder: {folder}',
    '[WARN] В файле нет SNP.': '[WARN] The file has no SNPs.',
    '[WARN] Ни один SNP не принят: проверьте геном и файл соответствий названий хромосом.':
        '[WARN] No SNP was accepted: check the genome and the chromosome name mapping file.',
    '[INFO] Сохранено: {out} · строк {rows} · SNP принято {accepted}, с кандидатами {candidates}, '
    'с robust-кандидатами {robust}':
        '[INFO] Saved: {out} · rows {rows} · SNPs accepted {accepted}, with candidates {candidates}, '
        'with robust candidates {robust}',
    # qtprimer3_vis, theme: the main window.
    'Системная': 'System',
    'Светлая': 'Light',
    'Тёмная': 'Dark',
    'Журнал расчёта': 'Run log',
    'Копировать': 'Copy',
    'Журнал': 'Log',
    'Тема: {mode} ▾': 'Theme: {mode} ▾',
    'Системный': 'System',
    'системный': 'system',
    'Язык: {name} ▾': 'Language: {name} ▾',
    'Язык': 'Language',
    'Язык изменится после перезапуска программы.': 'The language will change after the program is restarted.',
    'Новый расчёт': 'New run',
    'Ручной выбор у {count} SNP будет сброшен. Запустить?':
        ('The choice made by hand for {count} SNP will be reset. Run?',
         'Choices made by hand for {count} SNPs will be reset. Run?'),
    'Идёт расчёт…': 'Running…',
    'Останавливаю расчёт…': 'Stopping the run…',
    'Идёт расчёт… {done} из {count} SNP': ('Running… {done} of {count} SNP', 'Running… {done} of {count} SNPs'),
    'Ошибка расчёта: {message}': 'Run failed: {message}',
    'Расчёт остановлен': 'Run stopped',
    'Выберите файл SNP и геном слева и нажмите «Запустить»':
        'Choose an SNP file and a genome on the left and press “Run”',
    'В файле нет SNP': 'The file has no SNPs',
    'Готово: {count} SNP · {parts}': ('Done: {count} SNP · {parts}', 'Done: {count} SNPs · {parts}'),
    'Выгрузка': 'Export',
    'Сначала выполните расчёт.': 'Run the analysis first.',
    'Файл занят другой программой. Закройте его и повторите выгрузку.':
        'The file is open in another program. Close it and export again.',
    'Нет доступа к папке {folder}. Запустите программу из папки, куда разрешена запись.':
        'No access to the folder {folder}. Start the program from a folder where writing is allowed.',
    'Не удалось сохранить: {error}': 'Could not save: {error}',
    'Сохранено: {path}': 'Saved: {path}',
    'Результаты сохранены:\n{path}': 'Results saved:\n{path}',
    'Ошибка: {error}. Подробности — в журнале.': 'Error: {error}. Details are in the log.',
    'Ошибка': 'Error',
    '{error}\n\nПодробности — в окне «Журнал».': '{error}\n\nDetails are in the “Log” window.',
    # launch_panel: the «Запуск» column and the advanced dialog.
    'Нет': 'Off',
    'Если нет надёжной схемы': 'If there is no robust assay',
    'Для каждого SNP': 'For every SNP',
    'Выбрать файл…': 'Choose a file…',
    'Выберите файл SNP': 'Choose an SNP file',
    'Выберите геном': 'Choose a genome',
    'Файл не найден': 'File not found',
    'Включая непродаваемые': 'Including non-commercial',
    'Все, включая непродаваемые': 'All, including non-commercial',
    'Все продаваемые': 'All commercial',
    '{names} и ещё {count}': '{names} and {count} more',
    'Видимые фрагменты от, п.н.': 'Visible fragments from, bp',
    'Нижняя граница видимых полос: более короткие фрагменты на геле не видны и в модели скрыты.':
        'The lower limit of visible bands: shorter fragments are not seen on the gel and are hidden in the model.',
    'до, п.н.': 'to, bp',
    'Верхняя граница видимых полос: более длинные фрагменты в модели скрыты.':
        'The upper limit of visible bands: longer fragments are hidden in the model.',
    'Наименьшее различие видимых полос REF и ALT.': 'The smallest difference between visible REF and ALT bands.',
    'Разрешение, %': 'Resolution, %',
    'Модель полос без учёта яркости. Различие должно превышать максимум из порога в п.н. '
    'и процента длины. Параметры требуют проверки на вашем геле.':
        'A band model without brightness. A difference must exceed the larger of the bp threshold '
        'and the percentage of the length. Check the parameters on your own gel.',
    'Фланк, п.н.': 'Flank, bp',
    'Число оснований по обе стороны SNP.': 'Number of bases on each side of the SNP.',
    'Схем на SNP': 'Assays per SNP',
    '0 — все схемы.': '0: all assays.',
    'Макс. разрезов': 'Max. cuts',
    'Объединять изошизомеры': 'Group isoschizomers',
    'Одна строка на сайт узнавания для SNP: лучшая схема, остальные ферменты с тем же '
    'сайтом — в столбце «Изошизомеры».':
        'One row per recognition site for an SNP: the best assay, with the other enzymes of the same '
        'site in the “Isoschizomers” column.',
    'SNP от края ампликона, п.н.': 'SNP from amplicon end, bp',
    'Primer3 ставит SNP не ближе заданного числа п.н. к каждому концу ампликона, чтобы '
    'диагностический фрагмент был видим. Значение должно быть не меньше мин. фрагмента. '
    'Не действует в поиске по окну.':
        'Primer3 keeps the SNP at least this many bp from each amplicon end, so that the '
        'diagnostic fragment is visible. The value should be no less than the min. fragment. '
        'Has no effect in the window search.',
    'Исключать сайт фермента в праймере': 'Exclude an enzyme site in a primer',
    'Сайт внутри праймера есть в продукте при любом аллеле, а разрез у самого конца ДНК '
    'ненадёжен.':
        'A site inside a primer is in the product with either allele, and a cut at the very end of the DNA '
        'is unreliable.',
    'Учитывать повторы (soft-mask)': 'Account for repeats (soft-mask)',
    'Строчные буквы FASTA (повторы в геномах NCBI/Ensembl) не допускаются на 3′-конце праймера.':
        'Lower-case FASTA letters (repeats in NCBI/Ensembl genomes) are not allowed at the primer 3′ end.',
    'Исключать схемы с фоновыми вариантами': 'Exclude assays with background variants',
    'Схемы, где известные варианты из VCF меняют разрезы или длины полос.':
        'Assays where known variants from the VCF change cuts or band sizes.',
    'dCAPS (праймер с несовпадением)': 'dCAPS (mismatch primer)',
    'Несовпадений, до': 'Mismatches, up to',
    'Только для схем dCAPS: наибольшее число несовпадений в праймере.':
        'dCAPS assays only: the largest number of mismatches in the primer.',
    'Tm праймеров до, °C': 'Primer Tm up to, °C',
    'Только для схем dCAPS: верхняя граница Tm праймера с несовпадением и второго праймера '
    '(нижняя — Tm на панели «Запуск»). Длинный праймер отодвигает разрез от конца продукта '
    'и увеличивает различие полос.':
        'dCAPS assays only: the upper Tm limit of the mismatch primer and the second primer '
        '(the lower one is the Tm of the “Run” panel). A longer primer moves the cut away from the product end '
        'and enlarges the band difference.',
    'Продукт от, п.н.': 'Product from, bp',
    'Только для схем dCAPS: на коротком продукте различие ~20 п.н. хорошо видно на геле.':
        'dCAPS assays only: on a short product a ~20 bp difference is easy to see on a gel.',
    'Только для схем dCAPS: наибольший размер продукта.': 'dCAPS assays only: the largest product size.',
    'Только для схем dCAPS: наименьшее различие видимых полос REF и ALT; заменяет общий порог Δ.':
        'dCAPS assays only: the smallest difference between visible REF and ALT bands; '
        'replaces the general Δ threshold.',
    'Без Primer3 анализируется окно вокруг SNP.': 'Without Primer3 the window around the SNP is analysed.',
    'Только gain/loss': 'Gain/loss only',
    'Пропускать нормализацию': 'Skip normalization',
    'Отключает только перестановку REF/ALT.': 'Turns off only the REF/ALT swap.',
    'Не сохранять нормализацию': 'Do not save normalization',
    'Гель': 'Gel',
    'Поиск': 'Search',
    'Проверки': 'Checks',
    'Особые случаи': 'Special cases',
    'Дополнительные параметры': 'Advanced parameters',
    'Все': 'All',
    'Не ограничивать': 'No limit',
    'Сбросить к рекомендуемым': 'Reset to recommended',
    'Готово': 'Done',
    'Запуск': 'Run',
    'Свернуть панель': 'Collapse the panel',
    'Данные': 'Data',
    'CSV или TSV: chrom, pos, ref, alt': 'CSV or TSV: chrom, pos, ref, alt',
    'Геном': 'Genome',
    'изменить': 'change',
    'Без соответствий': 'No mapping',
    'Автоматически': 'Automatic',
    'не обязательно': 'optional',
    'VCF или VCF.gz известных вариантов породы или популяции.':
        'VCF or VCF.gz of known variants of the breed or population.',
    'Ферменты': 'Enzymes',
    'dCAPS: праймер с 1–2 несовпадениями у 3′-конца создаёт сайт фермента только '
    'в одном аллеле. Различие полос ≈ длине праймера.':
        'dCAPS: a primer with 1–2 mismatches near the 3′ end creates the enzyme site in one allele '
        'only. The band difference ≈ the primer length.',
    'Дополнительно…': 'Advanced…',
    'Запустить': 'Run',
    'Развернуть панель «Запуск»': 'Expand the “Run” panel',
    'Геном FASTA': 'FASTA genome',
    'FASTA (*.fa *.fasta *.fna);;Все файлы (*)': 'FASTA (*.fa *.fasta *.fna);;All files (*)',
    'Хромосомы: {name}': 'Chromosomes: {name}',
    'Хромосомы: как в файле SNP': 'Chromosomes: as in the SNP file',
    'Хромосомы: {name} (авто)': 'Chromosomes: {name} (auto)',
    'Файл соответствий хромосом': 'Chromosome name mapping file',
    'Текст (*.txt);;Все файлы (*)': 'Text (*.txt);;All files (*)',
    'Файл SNP': 'SNP file',
    'Таблицы (*.tsv *.csv *.txt);;Все файлы (*)': 'Tables (*.tsv *.csv *.txt);;All files (*)',
    'VCF известных вариантов': 'VCF of known variants',
    'VCF (*.vcf *.vcf.gz *.gz);;Все файлы (*)': 'VCF (*.vcf *.vcf.gz *.gz);;All files (*)',
    'ошибка': 'error',
    'изменено {count}': 'changed {count}',
    'Остановить': 'Stop',
    'Подготовка…': 'Preparing…',
    '{done} из {count} SNP': ('{done} of {count} SNP', '{done} of {count} SNPs'),
    # snp_list: the SNP list.
    'Полосы REF | ALT': 'Bands REF | ALT',
    'Продукт': 'Product',
    'Схем': 'Assays',
    'Все {count}': 'All {count}',
    'Выгрузить': 'Export',
    'Excel: все SNP': 'Excel: all SNPs',
    'Excel: только показанные ({count})': 'Excel: shown only ({count})',
    'CSV: все схемы': 'CSV: all assays',
    'Поиск по ID, координате или ферменту': 'Search by ID, position or enzyme',
    'Сбросить фильтр': 'Reset filter',
    'Столбцы…': 'Columns…',
    'как в файле': 'as in file',
    'по столбцу «{column}»': 'by column “{column}”',
    'Показано {shown} из {total} · порядок: {order}': 'Shown {shown} of {total} · order: {order}',
    'Нет SNP с такими условиями': 'No SNPs match these conditions',
    # details_panel: the details column and the primer card.
    'Выберите SNP в списке': 'Choose an SNP in the list',
    '+{count} изошизомер': '+{count} isoschizomer',
    '+{count} изошизомеров': '+{count} isoschizomers',
    'Карта': 'Map',
    '⤢ Крупно': '⤢ Enlarge',
    'Карта с последовательностями и гель в отдельном окне, с сохранением PNG и SVG':
        'The map with sequences and the gel in a separate window, with PNG and SVG export',
    'Копировать для заказа': 'Copy for ordering',
    'Две строки «имя ⇥ последовательность» для бланка заказа синтеза':
        'Two lines “name ⇥ sequence” for an oligo synthesis order form',
    'Копировать всё': 'Copy all',
    'Текстовая сводка схемы для лабораторного журнала': 'A text summary of the assay for the lab notebook',
    'строка {line} файла': 'line {line} of the file',
    'Схемы ({count}) · ★ идёт в Excel': 'Assays ({count}) · ★ goes to Excel',
    'Выбрана для Excel': 'Chosen for Excel',
    'Сделать выбранной': 'Make it the chosen one',
    'лучшая': 'best',
    '{ref_ref} | {alt_alt} · продукт {size}': '{ref_ref} | {alt_alt} · product {size}',
    'Изошизомеры: {names}': 'Isoschizomers: {names}',
    '{length} нт · Tm {tm} °C': '{length} nt · Tm {tm} °C',
    'Продукт {size} п.н. · {enzyme} ({site})': 'Product {size} bp · {enzyme} ({site})',
    'Поставщики: {suppliers}': 'Suppliers: {suppliers}',
    'Строчные буквы — намеренные несовпадения.': 'Lower-case letters are intended mismatches.',
    'Эта схема не выбрана: в Excel идёт ★ {enzyme}. Отметьте её звёздочкой, чтобы копировать.':
        'This assay is not chosen: ★ {enzyme} goes to Excel. Star it to copy.',
    # visualization: the amplicon map and the virtual gel.
    'Карта ампликона': 'Amplicon map',
    'Уменьшить схему': 'Zoom out',
    'Увеличить схему; Ctrl + колесо мыши': 'Zoom in; Ctrl + mouse wheel',
    'Вписать': 'Fit',
    'Сбросить масштаб и показать всю схему': 'Reset the zoom and show the whole map',
    'Сохранить PNG': 'Save PNG',
    'Вся схема в высоком разрешении': 'The whole map at high resolution',
    'Сохранить SVG': 'Save SVG',
    'Векторная схема для печати и редактирования': 'A vector map for printing and editing',
    'Последовательности 5′ → 3′': 'Sequences 5′ → 3′',
    'Копировать REF': 'Copy REF',
    'Полная последовательность REF без разметки': 'The full REF sequence without marks',
    'Копировать ALT': 'Copy ALT',
    'Полная последовательность ALT без разметки': 'The full ALT sequence without marks',
    'Копировать праймеры': 'Copy primers',
    'Оба праймера в направлении 5′ → 3′': 'Both primers in the 5′ → 3′ direction',
    'Ctrl + колесо — масштаб; перетаскивание — перемещение схемы.': 'Ctrl + wheel: zoom; drag: move the map.',
    'Выберите результат в таблице, чтобы увидеть ампликон и сравнить REF с ALT.':
        'Choose a result to see the amplicon and compare REF with ALT.',
    'Карта появится после выбора результата': 'The map appears when a result is chosen',
    'ALT недоступна: позиция SNP отсутствует или находится за пределами ампликона.':
        'ALT is not available: the SNP position is missing or outside the amplicon.',
    'ALT недоступна: для двух дорожек требуется одиночная замена A/C/G/T.':
        'ALT is not available: two lanes need a single A/C/G/T substitution.',
    'ALT недоступна: REF {ref} не совпадает с основанием {base} в FASTA.':
        'ALT is not available: REF {ref} does not match the base {base} in the FASTA.',
    'Не удалось получить последовательность ампликона.': 'Could not get the amplicon sequence.',
    'Левый праймер не найден в ампликоне.': 'The left primer is not found in the amplicon.',
    'Сайт посадки правого праймера не найден в ампликоне.':
        'The binding site of the right primer is not found in the amplicon.',
    'Фермент не найден в базе Biopython; разрезы не показаны.':
        'The enzyme is not in the Biopython database; cuts are not shown.',
    'Длина последовательности не совпадает с геномными координатами.':
        'The sequence length does not match the genome coordinates.',
    'Сравните разрезы и длины фрагментов; наведите курсор на метку для точных координат.':
        'Compare the cuts and fragment sizes; hover over a mark for exact coordinates.',
    'Фермент не выбран': 'No enzyme chosen',
    '{chrom}:{start}–{end}  ·  {length} п.н.  ·  геномные координаты 1-based':
        '{chrom}:{start}–{end}  ·  {length} bp  ·  1-based genome coordinates',
    'Появление / исчезновение сайта': 'Site gain / loss',
    'Изменение фрагментов': 'Fragment change',
    '{label}  5′ → 3′  ·  {length} п.н.': '{label}  5′ → 3′  ·  {length} bp',
    '{label}: не подобран': '{label}: not designed',
    '{label}: {primer}\nПосадка: {start}–{end} в ампликоне (1-based)':
        '{label}: {primer}\nBinds at {start}–{end} in the amplicon (1-based)',
    'Левый праймер →': 'Left primer →',
    '← Правый праймер': '← Right primer',
    'Дорожка недоступна — см. сообщение над схемой': 'The lane is not available; see the message above the map',
    '{label}: после {position} п.н. от начала ампликона': '{label}: after {position} bp from the amplicon start',
    'Граница ампликона: {position} п.н.': 'Amplicon end: {position} bp',
    'Фрагменты: {detail} п.н.': 'Fragments: {detail} bp',
    'Разрезы не рассчитаны: фермент не выбран или неизвестен.':
        'Cuts not calculated: no enzyme chosen or the enzyme is unknown.',
    'Фрагментов: {count}; полные длины — в последовательностях и подсказке.':
        'Fragments: {count}; the full sizes are in the sequences and the tooltip.',
    'SNP: основание {base} в ампликоне; геномная позиция {position}':
        'SNP: base {base} of the amplicon; genome position {position}',
    'Шкала и метки разрезов: п.н. от начала ампликона · длины фрагментов указаны под дорожками':
        'Scale and cut marks: bp from the amplicon start · fragment sizes are under the lanes',
    'SNP  {ref} → {alt}  ·  основание {base}': 'SNP  {ref} → {alt}  ·  base {base}',
    'Основания {start}–{stop} ампликона (1-based)': 'Bases {start}–{stop} of the amplicon (1-based)',
    'REF: референсная последовательность · ALT: последовательность после замены':
        'REF: the reference sequence · ALT: the sequence after the substitution',
    'Показаны основания {begin}–{end} из {length}. Кнопки копируют полную последовательность.':
        'Bases {begin}–{end} of {length} are shown. The buttons copy the full sequence.',
    'Левый  5′→3′: {left}\nПравый 5′→3′: {right}': 'Left   5′→3′: {left}\nRight  5′→3′: {right}',
    '{label}, фрагменты: {fragments} п.н.': '{label}, fragments: {fragments} bp',
    'Ампликон {start}–{stop} · геном {genome_start}–{genome_stop}':
        'Amplicon {start}–{stop} · genome {genome_start}–{genome_stop}',
    '{allele}: скопировано {length} п.н.': '{allele}: {length} bp copied',
    'Праймеры скопированы в формате FASTA (5′ → 3′).': 'Primers copied in FASTA format (5′ → 3′).',
    'Сначала выберите результат для визуализации.': 'Choose a result to draw first.',
    'Ампликон, праймеры 5′→3′, SNP и разрезы рестриктазы. Координаты генома 1-based; '
    'разрезы — п.н. от начала ампликона.':
        'Amplicon, primers 5′→3′, the SNP and restriction cuts. Genome coordinates are 1-based; '
        'cuts are bp from the amplicon start.',
    'Сохранить карту ампликона': 'Save the amplicon map',
    'Экспорт схемы': 'Map export',
    'Без фермента': 'No enzyme',
    'Масштаб экспорта должен быть от 0.5 до 8.': 'The export scale must be between 0.5 and 8.',
    'Не удалось сохранить PNG: {path}': 'Could not save the PNG: {path}',
    'Не удалось сохранить SVG: {error}': 'Could not save the SVG: {error}',
    'Не удалось начать экспорт SVG: {path}': 'Could not start the SVG export: {path}',
    'Показать весь гель': 'Show the whole gel',
    'Гель в высоком разрешении': 'The gel at high resolution',
    'Векторный рисунок геля': 'A vector picture of the gel',
    'Выберите результат в таблице, чтобы увидеть ожидаемые полосы.': 'Choose a result to see the expected bands.',
    'Длина продукта неизвестна.': 'The product size is unknown.',
    'Виртуальный гель': 'Virtual gel',
    'Продукт {product} п.н. · видимые фрагменты {low}–{high} п.н. · полная рестрикция':
        'Product {product} bp · visible fragments {low}–{high} bp · complete digestion',
    'Маркер: {length} п.н.': 'Marker: {length} bp',
    '{label}: {length} п.н., доля ДНК {share}': '{label}: {length} bp, DNA share {share}',
    ' — вне видимого диапазона': ', outside the visible range',
    '{fragments} п.н.': '{fragments} bp',
    'Схема по модели программы: яркость пропорциональна массе фрагмента, у гетерозиготы — '
    'по половине каждого аллеля; в скобках и бледно — фрагменты вне видимого диапазона.':
        'Drawn from the program model: brightness follows fragment mass, a heterozygote has half of each '
        'allele; fragments outside the visible range are faint and in brackets.',
    'Виртуальный гель RFLP': 'RFLP virtual gel',
    'Ожидаемые полосы REF/REF, REF/ALT и ALT/ALT по модели RFLP Picker.':
        'Expected REF/REF, REF/ALT and ALT/ALT bands by the RFLP Picker model.',
    'Сохранить гель': 'Save the gel',
    'Экспорт геля': 'Gel export',
    'Карта ампликона · RFLP REF / ALT': 'Amplicon map · RFLP REF / ALT',
    'Закрыть': 'Close',
    'Гель недоступен: нет последовательности ALT (см. карту ампликона).':
        'The gel is not available: there is no ALT sequence (see the amplicon map).',
    'Не задан FASTA для визуализации. Запустите анализ заново.': 'No FASTA set for drawing. Run the analysis again.',
    'В выбранной строке нет координат ампликона.': 'The chosen row has no amplicon coordinates.',
    'Некорректные координаты ампликона: {region}': 'Invalid amplicon coordinates: {region}',
    'Не удалось открыть FASTA: {path}\n{error}': 'Could not open the FASTA: {path}\n{error}',
    'Не удалось получить последовательность из FASTA: {region}\n{error}':
        'Could not read the sequence from the FASTA: {region}\n{error}',
}

CORE = [
    # Reasons in the result rows.
    (r'SNP ближе (\d+) п\.н\. к краю окна; увеличьте фланг',
     r'SNP is closer than \1 bp to the window edge; increase the flank'),
    (r'Primer3 не вернул подходящих пар', 'Primer3 returned no suitable pairs'),
    (r'Не найдено сочетание для различения трёх генотипов',
     'No primer–enzyme combination distinguishes the three genotypes'),
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
    (r'Primer3 вернул праймеры вне окна или перекрывающие SNP',
     'Primer3 returned primers outside the window or over the SNP'),
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
    (r'Минимальный размер продукта dCAPS больше максимального', 'Minimum dCAPS product size is greater than maximum'),
    (r'Минимальный Tm праймеров dCAPS больше максимального', 'Minimum dCAPS primer Tm is greater than maximum'),
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
    (r'Оценка полос: модель присутствия, полная рестрикция; яркость не оценивается\.',
     'Band model: presence of bands, complete digestion; band brightness is not assessed.'),
    (r'пары праймеров с посторонним продуктом ПЦР исключаются', 'primer pairs with another PCR product are excluded'),
    (r'пары праймеров с посторонним продуктом ПЦР отмечаются', 'primer pairs with another PCR product are flagged'),
    (r'\(до (\d+) несовпадений вне 3′-конца, продукт до (\d+) п\.н\.\)',
     r'(up to \1 mismatches outside the 3′ end, product up to \2 bp)'),
    (r'Специфичность: проход по геному для (\d+) праймеров', r'Specificity: a pass over the genome for \1 primers'),
    (r'Специфичность: просмотрено (\d+)% генома', r'Specificity: \1% of the genome scanned'),
    (r'Все пары праймеров дают посторонние продукты ПЦР или лежат в повторах генома',
     'All primer pairs give other PCR products or lie in genome repeats'),
    (r'\[INFO\] Отбор: ', '[INFO] Filters: '),
    (r'только продаваемые ферменты', 'commercial enzymes only'),
    (r'изошизомеры объединены \(одна строка на сайт\)', 'isoschizomers grouped (one row per site)'),
    (r'без сайта фермента в праймерах', 'no enzyme site in primers'),
    (r'повторы \(строчные буквы FASTA\) не на 3′-концах праймеров', 'repeats (lower-case FASTA) not at primer 3′ ends'),
    (r'без фоновых вариантов, меняющих разрезы', 'no background variants that change cuts'),
    (r'SNP не ближе (\d+) п\.н\. к концам ампликона', r'SNP at least \1 bp from amplicon ends'),
    (r'dCAPS если нет robust-схемы, несовпадений до (\d+)', r'dCAPS if no robust assay, mismatches up to \1'),
    (r'dCAPS для каждого SNP, несовпадений до (\d+)', r'dCAPS for every SNP, mismatches up to \1'),
    (r', Tm до ([\d.]+) °C, продукт (\d+)–(\d+) п\.н\., Δ (\d+) п\.н\.',
     r', Tm up to \1 °C, product \2–\3 bp, Δ \4 bp'),
    (r'Фоновые варианты: (\d+) в окнах SNP из (\d+) записей VCF',
     r'Background variants: \1 in SNP windows out of \2 VCF records'),
    (r': пар (\d+), кандидатов (\d+), неразличимых (\d+)', r': pairs \1, candidates \2, indistinguishable \3'),
    (r'В окнах SNP нет строчных букв: FASTA без soft-mask, маскирование повторов не действует\.',
     'No lower-case letters in SNP windows: the FASTA is not soft-masked, repeat masking has no effect.'),
    (r'NORMALIZE: загружено вариантов: ', 'NORMALIZE: variants loaded: '),
    # Logger lines; the command line shows them.
    (r'Файл names не найден: ', 'Names file not found: '),
    (r'Ошибка при загрузке names: ', 'Error loading the names file: '),
    (r'Файл с вариантами не найден: ', 'Variants file not found: '),
    (r'Ошибка парсинга строки: ', 'Line parse error: '),
    (r'Ошибка при чтении файла с вариантами: ', 'Error reading the variants file: '),
]
