# Mat_belle

[Ссылка на структуру данных](https://belle.kek.jp/~hastings/secured/pantherdoc/b20030807_1600/)
## Первый эксперимент

Сбор данных:

* Оставил только полулептонные распады lamc_pm

* Добавил распады для D0, D_pm

* Изменил заполнение на каналов на map

* Убрал каты лептонов

Результат:

* В тагитированной массе lamc видно порго масс нейтрино

* Решили обрезать нейтрино m2_nu  (-2, 4)

## Эксперимент второй

Сбор данных:

* Добавил эталонный канал l -> ap, k_p, pi_m (чистый)

* Добавил каналы подобные полулептоным (каналы 3, 4 в lamct)

* Добавил отбор событий (посже исправил and -> or)

Результат:

* Отбор событий с катими $((p_{event} - (p_u - p_{\Lambda}))^2 > 3.5^2 \lor  (N_{ntr} \ge 1))$ обрезает статистику на ~ 99.99%

* **В сборе была ошибка канал 4 не собирался $\implies$ процент меньше**

## Эксперимент третий

Задачи:

* [x] Добавить в события Ups каналы с rho4

* [ ] Добавить [$D \to \bar K K \bar K$](https://pdg.lbl.gov/2020/tables/rpp2020-tab-mesons-charm.pdf) (канал имет вероятность $\sim 10^{-3} \% $, что мало)

* Добавить каналы $D_s^\pm$  
    * [x] $D_s^\pm \to K^+ K^- \pi^\pm$
    * [x] $D_s^\pm \to K_S \pi^\pm$

* Добавить каналы U
    * [x] $U \to \Lambda_c^\mp D_s^\pm \bar K^\mp p/\bar p$
    * [x] $U \to \Lambda_c^\pm D_s^\mp \bar \Lambda / \Lambda$
    * [x] $U \to \Lambda_c^\mp D^\pm (p/\bar p) K^\mp (K^0/\bar K^0) + \pi^0$ нужно собрать $K^0$
    * [x] $U \to \Lambda_c^\pm (\bar D^0/D^0)  K^\mp (\bar \Lambda/\Lambda)  + \pi^0$

* [x] Добавить переменную сум. en фотонов 

* [x] Записать импульс лептона


## Эксперимент 4

* [x] Переделать эенергию фотонов checkSame
* [x] Добавить rho6, не делать тк много комбинаций
* [ ] Поробовать заменить $\pi^+ \pi^- $ на $K^+ K^-$ и  $\pi^\pm \pi^0 $ на $K^\pm K_S$ 
* [x] Добавить каналы $D^\pm \to K^\mp 2\pi^\pm \pi^0, \ D^\pm \to K_S \pi^\pm \pi^0$ 
* Еще поменял файл test.cc, теперь *args к работает как dictionary.
* Добавил автоматическую подгонку границ гистораммы.
* Добавил вывод нескольких графиков в один файл.
* Добавил фиты, адаптировал их.
* Переписал Makefile теперь можно работаь с любым С файлом.

## Эксперимент 5

* Убрал все собыбытия с $pi^0$
* Добавил $D^*$ и каналы $\Lambda p D^{*+} \pi^-, \Lambda p D^{*0}$
* Убрал все собыбытия, которые могли вдть двойноый счет с события ми из $D^{*0}$

## На данный момент

$e^+ e^- \to \Lambda_c^- X_c$
* [1] $X_c \to D^0 p$
* [2] $X_c \to D^+ p  pi^-$
* [3] $X_c \to D^{*+} p  pi^-$
* [4] $X_c \to D^{*0} p$

./cut "out = sqrt(rm2l)" "cut = abs(mach) < 0.01 && abs(ml - 2.28646) < 0.015 && chl == 5 && abs(rm2n) < 1" "chu = 1-7" "down = 1" "up = 4" "nbins = 50"

./cut "out = sqrt(rm2l)" "cut =   abs(mach) < 0.01 && abs(ml) < 0.015 && chl == 5 && abs(rm2n) < 1 && (abs(dm_dst - 0.142014) < 0.003 || chxc != 4)" "chu = chxc = 1-4" "down = 1" "up = 4" "nbins = 50" "fname = rm_lam_all_dt"

./cut "out = sqrt(rm2l)" "cut = chl == 5 && abs(rm2n) < 0.01 && abs(ml) < 0.015 && (abs(mach) < 0.015 || chxc <= 2) &&  (abs(mach) < 0.03 || chxc >= 3) " "chu = chxc = 1-4" "down = 1.5" "up = 3.5" "nbins = 50" "fname = rm_lam_all_dt_new"

./cut "out = 3.141592 - acos(cos_lam_)" "cut = rm2l > 0 &&  abs(sqrt(rm2l) - 2.28646) < 0.1 && chl == 5 && abs(rm2n) < 0.01 && abs(ml) < 0.015 && (abs(mach) < 0.015 || chxc <= 2) &&  (abs(mach) < 0.03 || chxc >= 3) " "chu = chxc = 1-4" "down = 0.0000000000001" "up = 0.01" "nbins = 50" "fname = p_angl_tag_taging"