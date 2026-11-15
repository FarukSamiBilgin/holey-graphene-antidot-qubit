> Yazarın Word raporundan (`holey_graphene_qubit_negative_result.docx`) Markdown'a dönüştürülmüştür. Converted from the author's Turkish negative-result report. An English summary is in the repository README.

**Delikli-grafen kusur-qubitinde parazit kenar momenti ile qubit kalitesi arasındaki ödünleşme ölçek-kararlı değildir: önceden-kayıtlı bir negatif sonuç**

*A pre-registered negative result on the parasitic-edge-moment versus qubit-quality trade-off in holey-graphene defect qubits*

F. Sami

Bağımsız hesaplamalı çalışma · 2026-05-19

**Özet**

Delikli grafen (graphene antidot lattice) içinde, antidot dizisinden bir deliğin kasıtlı olarak atlanmasıyla tanımlanan kusur durumu, uzun süredir katı-hal tabanlı bir spin-qubit adayı olarak önerilmektedir (Pedersen ve ark., 2008). Ayrı bir çalışma kümesi (Vanević ve ark., 2009) ise antidot kenarlarındaki yerel manyetik momentlerin, tasarlanan kusur durumuyla yarışan güçlü bir dekoherans kaynağı olduğunu işaret etmiştir. Bu iki gözlem literatürde hiç aynı geometri taramasında karşı karşıya getirilmemiştir: bir geometri ailesi boyunca qubit kalitesi ile parazit moment ÖRTÜŞÜYOR mu, yoksa ÖDÜNLEŞİYOR mu? Bu çalışma tam olarak bu soruyu, periyodik ölçülebilir-uyumlu bir sıkı-bağlanım + ortalama-alan Hubbard hattı kurarak, sonuç görülmeden sabitlenmiş bir karar ölçütüyle test eder. Küçük süper-hücrede (Nd=3) ve Γ-noktasında görülen bir ödünleşme eğilimi, önceden-kayıtlı dayanıklılık testi uygulandığında — k-örnekleme ve süper-hücre boyutunun büyütülmesi (Nd=4) altında — korunmamaktadır: işaret sırası tersine döner. Önceden tanımlanmış ölçüte göre sonuç bir ARTIFACT'tir: bu yaklaşımda bu ölçekte görülen ödünleşme, sonlu-hücre / Γ-noktası seçiminin bir yan ürünüdür, fiziksel bir etki değildir. Bu, Vanević sorusuna verilmiş dürüst, nicel ve sınırlayıcı bir cevaptır ve bir metodolojik not olarak raporlanmaktadır.

## 1. Giriş ve problem tanımı

Grafen antidot örgüleri (delikli grafen), grafen tabakasından periyodik bir delik dizisi oyularak elde edilen ve bu sayede ayarlanabilir bir bant aralığı kazanan bir malzeme sınıfıdır. Pedersen ve arkadaşları (2008), bu dizide tek bir deliğin oyulmadan bırakılmasının, aralık içinde uzaysal olarak yerelleşmiş ayrık bir elektronik durum yarattığını ve bunun bir spin-qubit taşıyıcısı olarak kullanılabileceğini öne sürmüştür. Sonraki bir çalışma (Pedersen grubu, 2013) periyodik örgünün spin-polarizasyonunu ve aralığını geometriye göre sistematik olarak taramış; ancak bu, kusur-qubitin kendisini değil, kusursuz periyodik örgüyü incelemiştir.

Bağımsız olarak, Vanević ve arkadaşları (2009), antidot kenarlarındaki ve aralık-ortası (midgap) durumlardaki yerel manyetik momentlerin, tasarlanan kusur durumuyla yarışan ve onu kirletebilen güçlü bir defaz (dephasing) kaynağı olduğunu teşhis etmiştir — fakat bu çalışma bir çözüm ya da geometri taraması sunmamıştır.

Bu çalışmanın hedeflediği özgün boşluk şudur: **Pedersen-2008'in kusur-qubit NESNESİNİ, Pedersen-2013'ün spin-polarize tarama YÖNTEMİNİ ve Vanević-2009'un yarışma SORUSUNU tek bir taramada birleştirmek**. Yani her geometri için qubit kalitesini ve parazit momenti AYNI hesapta ölçmek ve şu soruyu sormak: bu iki büyüklük örtüşür mü (literatürün ‘armchair kenar = momentsiz, bedava’ anlatısı doğru) yoksa ödünleşir mi (anlatıyı sorgulayan, Vanević yarışmasını nicelleştiren bir sonuç)?

Önemli bir metodolojik nokta baştan kabul edilmiştir: ‘armchair kenar moment üretmez’ kuralı iyi bilindiğinden, temiz bir pozitif sonuç özgün olmazdı. Bilimsel değeri olan çıktı, bir ödünleşmenin ölçek değişimi altında HAYATTA KALIP kalmadığını haritalamaktır — ve dayanıklı bir negatif de en az pozitif kadar yayımlanabilir bir sonuçtur. Bu çalışma, bu disiplini sonuna kadar uygulamakta ve negatif sonucu ana bulgu olarak raporlamaktadır.

## 2. Literatürdeki konum ve özgünlük

İlgili çalışmalar ve bu çalışmanın onlardan ayrıldığı nokta:

- Pedersen ve ark., PRL 100, 136804 (2008): kurucu öneri; kusur-qubiti antidot atlanarak tanımlar. Saf sıkı-bağlanım, spin polarizasyonu yok.

- Pedersen grubu, PRB 88, 195418 (2013): DFT + ortalama-alan Hubbard ile spin-polarizasyon/aralık geometri taraması — fakat PERİYODİK örgünün, kusur-qubitin değil.

- Vanević ve ark. (2009): kenar/midgap yerel momentlerinin kusur durumuyla yarışan bir dekoherans kaynağı olduğunu teşhis eder — çözüm ya da tarama önermez.

- Twamley hattı (Cirio/Brennen/Twamley 2012; Johnsson/Brennen/Twamley 2016): levite süperiletken mikroküreler; grafen yalnızca mekanik bir öğe. (Bu çalışmanın başında levitasyon kapsam dışı bırakılmıştır: levitasyon mekanik/kenetlenme dekoheransını kaldırır, fakat asıl sorun olan elektronik kenar-moment defazını kaldırmaz.)

**Boşluk:** Hiçbir mevcut çalışma, kusur-qubitin kalitesi ile parazit kenar momentini aynı geometri taramasında, üstelik ölçek-dayanıklılığı önceden-kayıtlı bir ölçütle test edilmiş halde karşılaştırmamıştır. Bu çalışmanın katkısı tam olarak budur.

## 3. Yöntem

### 3.1 Geometri: periyodik ölçülebilir-uyumlu süper-hücre

Sonlu pul (flake) yaklaşımı, dikdörtgen pulun kendi sınırının kusur durumunu bastırması nedeniyle terk edilmiştir (bkz. §4, L1). Onun yerine, Pedersen 2008/2013 kurgusu izlenmiştir: grafen Bravais vektörleri ile bir antidot hücresi L×L grafen hücresi büyüklüğünde tanımlanır; kusur süper-hücresi Nd×Nd antidot hücresinden oluşur ve merkezi delik OYULMADAN bırakılır. Tüm yapı tam-sayı ölçülebilir-uyumludur, dolayısıyla periyodik sınır koşulları tamdır. C–C bağ uzunluğu birim alınmış, en-yakın-komşu hopping büyüklüğü t = 2.7 eV kullanılmıştır.

### 3.2 Qubit ekseni: sıkı-bağlanım, Fermi-civarı spektrum

Bloch Hamiltoniyeni seyrek matris olarak kurulur; sıfır enerji civarındaki durumlar shift-invert ARPACK (eigsh) ile çözülür. Bipartit yapıda tam sıfır modları bulunduğundan, σ=0 tekil olur; 1e-3 eV'lik küçük bir kaydırma (her antidot aralığının çok içinde) koşullanmayı düzeltir. Determinizm için sabit bir başlangıç vektörü kullanılır ve bant aralığı yalnızca dispersif bantlar arasında ölçülür; sıfır-mod sayısı ayrı raporlanır (bkz. §4, L3).

**Qubit kalite skoru: q = izolasyon × w_def**, burada izolasyon = aralık − |E_kusur| (kusur seviyesinin aralık içinde, bant kenarlarından ne kadar derinde olduğu) ve w_def = o seviyenin tasarlanan-kusur bölgesindeki ağırlığı. q maksimize edilmek istenir.

### 3.3 Parazit eksen: ortalama-alan Hubbard

Saf sıkı-bağlanım düzeyinde zigzag sıfır modları yalnızca dejenere durumlardır — manyetik moment DEĞİL. Moment (dolayısıyla projenin konusu olan parazit kenar momenti) ancak yerinde Coulomb U dahil edildiğinde ortaya çıkar. Collinear π-orbital ortalama-alan Hubbard, yarı-dolulukta öz-tutarlı çözülür. Yerel moment m_i = ⟨n_i↑⟩ − ⟨n_i↓⟩. Metrikler: μ_edge = jant atomu başına ortalama |m_i| (Vanević'in parazit kenar momenti), μ_def = kusur bölgesi atomu başına ortalama |m_i|, M_net = |Σ m_i| (küresel net spin; Lieb alt sınırı = |dN|). Σ|m_i| metriğinin neden yanlış olduğu §4'te (L4) açıklanmıştır.

**Parazit skoru: p = μ_edge** (minimize edilmek istenir). Ortalama-alan Hubbard U'ya duyarlı olduğundan tek bir kırılgan değer yerine bir U/t bandı kullanılır.

### 3.4 Önceden-kayıtlı dayanıklılık ölçütü

Bu çalışmanın bilimsel bütünlüğünün çekirdeği budur. Çekirdek tarama, Γ-noktasında ve tek bir Nd'de görülen tek bir korelasyona dayanır: alt-örgü dengesizliği |dN| büyüdükçe hem qubit izolasyonu q hem de parazit moment p birlikte büyür. Sorulacak tek soru, VERİ GÖRÜLMEDEN şöyle sabitlenmiştir:

*(|dN| ↑  ⇒  q ↑  VE  p ↑) sıralaması, (a) Γ yerine bir k-örneklemesi ve (b) süper-hücre boyutu Nd değişimi altında korunuyor mu?*

Karar kuralı (veriye göre ayarlanamaz biçimde önceden sabit):

- ROBUST: (q_yüksek − q_düşük) ve (p_yüksek − p_düşük) işaretinin TEST EDİLEN HER KOŞULDA (Γ ve k-örneklemesi, her Nd) pozitif kalması.

- ARTIFACT: işaretin HERHANGİ bir koşulda dönmesi.

Ek bağımsız dürüstlük bayrağı: kusur seviyesi k boyunca dağılırsa (E_disp > 0.10 × aralık), o geometride ‘qubit’ aslında bir banttır — yani q ekseninin kendisi o geometride bir Γ artefaktıdır ve açıkça raporlanır.

### 3.5 Hesaplama ortamı

Tüm hesap, yalnızca numpy + scipy + joblib ile, GPU olmadan, tek bir Python dosyasında (saf, paralelleştirilebilir fonksiyonlar) yürütülmüştür. Son dayanıklılık testi bir Intel i7-12700H / 16 GB RAM masaüstünde, kuantum_lab Python 3.9.25 ortamında (numpy 1.26.4, scipy 1.13.1, joblib 1.5.3) çalıştırılmıştır. Nd=4'te süper-hücre ~1700–2000 atoma ulaşmaktadır; bu, bu yaklaşımın masaüstü-uygulanabilir sınırına yakındır (bkz. §6).

## 4. Sekiz-katmanlı doğrulama zinciri ve yakalanan hatalar

Bu çalışmanın metodolojik gücü, her katmanın en küçük kesin sistemde, önceden tanımlanmış geçti/kaldı ölçütüyle doğrulanmasından ve — kritik biçimde — her katmanda yakalanan hatanın gizlenmeden raporlanmasından gelir. Hatasız sunulan bir hesaplama hattı değil, hatası bulunup düzeltilmiş ve izi şeffaf bırakılmış bir hat, çok daha güvenilirdir.

| Katman | Eklenen | Yakalanan hata | Düzeltme |
|---|---|---|---|
| L1 | Sonlu pul (flake) sıkı-bağlanım iskeleti | Dikdörtgen pulun kendi sınırı kusur durumunu bastırdı; kusur seviyesi pul kenarından ayırt edilemedi | Periyodik, ölçülebilir-uyumlu süper-hücreye geçildi (Pedersen kurgusu) |
| L2 | Periyodik ölçülebilir süper-hücre (graphene Bravais, antidot L×L, kusur Nd×Nd, merkezi delik atlanmış) | —  (doğrulama geçti: temiz ~1.5 eV TB aralığı; kusur durumu w_def≈0.36 ile çözüldü) | Katman onaylandı; üzerine inşa edildi |
| L3 | Determinizm: shift-invert eigsh sabit başlangıç vektörü; aralık dispersif bantlar arasında ölçülüyor | Tam sıfır modları nedeniyle çıktı tohuma (seed) bağlıydı; aynı girdi farklı sonuç veriyordu | Sabit v0 + sıfır-mod sayısı ayrı raporlanır; 3 koşumda bit-aynı kanıtlandı |
| L4 | Ortalama-alan Hubbard (collinear, yarı-dolu, Γ) | İlk metrik Σ\|m_i\| yanlıştı: kademeli AFM dalgalanmalarını topluyordu, momentsiz hücre bile ~2 puan alıyordu | μ_edge / μ_def (atom başına ortalama \|m\|) ve M_net = \|Σm_i\| metriklerine geçildi |
| L5 | Eşit-alan şekil normalizasyonu | Küçük R'de farklı şekiller aynı atom kümesine çöküyordu (sahte özdeşlik) | Her şekil aynı referans alanına ölçeklendi |
| L6 | Gapsiz geometrilerin ele alınışı | Gapsiz zigzag nan q veriyor ve sessizce eleniyordu; bu, verdict'i armchair lehine yanlı kılıyordu | nan q → q=0'a zorlandı (gerçekten en kötü qubit; karşılaştırmada tutuldu) |
| L7 | Sürekli kiralite χ denemesi → yapısal parmak-izi ile ayıklama | Ayrık örgüde küçük dönüşler aynı atomları siliyordu (χ donması); 20 'farklı' geometri aslında ~6 idi | Yapısal parmak-izi ile yinelenenler ayıklandı; sadece gerçekten farklı geometriler değerlendirildi |
| L8 | Dayanıklılık testi: k-örnekleme + çok-Nd, ÖNCEDEN-KAYITLI ölçüt | Bu katman testin kendisidir; ölçüt veri görülmeden sabitlendi | Sonuç: ARTIFACT (aşağıda) |

**Özellikle L7 (χ donması) doğrudan bu sonucu mümkün kıldı:** ayrık bir örgüde sürekli bir kiralite ekseni fiziksel olarak yoktur; küçük dönüşler aynı atomları siler. Bu fark edilmeseydi, ‘20 farklı geometride ödünleşme’ gibi sahte-güçlü ama aslında ~6 bağımsız geometriye dayanan bir iddia üretilebilirdi. Yapısal parmak-izi ayıklaması bunu engelledi ve taramayı dürüst kıldı.

## 5. Sonuçlar: önceden-kayıtlı dayanıklılık testi

Otomatik seçilen iki yapısal-olarak-farklı temsilci geometri: düşük dengesizlik LO = altıgen/armchair, R=2.8 (|dN|≈8) ve yüksek dengesizlik HI = altıgen/zigzag, R=2.8 (|dN|≈16). Her ikisi de Γ (k=1) ve 2×2 k-örneklemesinde, Nd=3 ve Nd=4'te, U/t = 1.15'te değerlendirilmiştir. Tam sayısal sonuçlar aşağıdadır (mavi satırlar Nd=3, kahverengi satırlar Nd=4):

| Nd | k | geometri | \|dN\| | gap | izol. | E_disp | w_def | μ_edge | q | p | düz? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | Γ | hex/armchair | 8 | 0.948 | 0.208 | 0.000 | 0.358 | 0.0298 | 0.0746 | 0.0298 | yes |
| 3 | Γ | hex/zigzag | 16 | 1.215 | 0.607 | 0.000 | 0.253 | 0.0962 | 0.1534 | 0.0962 | yes |
| 3 | 2×2 | hex/armchair | 8 | 0.941 | 0.187 | 1.466 | 0.244 | 0.0298 | 0.0457 | 0.0298 | NO |
| 3 | 2×2 | hex/zigzag | 16 | 1.211 | 0.590 | 0.000 | 0.277 | 0.0962 | 0.1633 | 0.0962 | yes |
| 4 | Γ | hex/armchair | 15 | 0.960 | 0.480 | 0.000 | 0.190 | 0.0307 | 0.0912 | 0.0307 | yes |
| 4 | Γ | hex/zigzag | 30 | nan | nan | 0.000 | 0.034 | 0.0967 | 0.0000 | 0.0967 | NO |
| 4 | 2×2 | hex/armchair | 15 | 0.964 | 0.482 | 0.964 | 0.196 | 0.0307 | 0.0946 | 0.0307 | NO |
| 4 | 2×2 | hex/zigzag | 30 | nan | nan | 0.000 | 0.030 | 0.0967 | 0.0000 | 0.0967 | NO |

**Nd=3 (her iki k koşulu): sıralama KORUNUR.** dq(yüksek−düşük) = +0.079 (Γ), +0.118 (2×2); dp = +0.066 (her ikisi). |dN| büyüdükçe hem q hem p büyür — çekirdek taramanın ödünleşme sinyali bu ölçekte mevcuttur.

**Nd=4 (her iki k koşulu): sıralama BOZULUR.** dq = −0.091 (Γ), −0.095 (2×2); dp = +0.066. q farkının işareti tersine döner. Dahası, Nd=4'te zigzag geometri tamamen gapsizdir (gap = nan, q = 0; 24 sıfır modu) ve armchair tarafında kusur seviyesi k altında güçlü dağılım gösterir (E_disp = 0.96, ‘düz? NO’) — yani o geometride qubit, lokalize bir durum değil bir banttır.

Bağımsız dağılım bayrağı dört koşulda tetiklenmiştir: (Nd=3, 2×2, armchair), (Nd=4, Γ, zigzag), (Nd=4, 2×2, armchair), (Nd=4, 2×2, zigzag). Yani q ekseninin Γ-noktası okuması, test edilen ölçeklerde kararlı bir lokalize-qubit yorumunu desteklememektedir.

**ÖNCEDEN-KAYITLI VERDICT: ARTIFACT. Sıralama Nd=4'e geçişte (hem Γ hem 2×2 k-örneklemesinde) tersine döner. Çekirdek taramanın ‘ödünleşme’ sinyali, sonlu-hücre / Γ-noktası seçiminin bir yan ürünüdür. Bu yöntem, bu ölçekte, etkiyi TESPİT EDEMEMEKTEDİR.**

*Bu verdict, sonuç görülmeden §3.4'te sabitlenen ölçütün doğrudan uygulanmasıdır; sonradan seçilmiş ya da veriye uydurulmuş bir yorum değildir.*

## 6. Tartışma ve dürüst sınırlar

Sonucun ne olduğu ve ne OLMADIĞI net ayrılmalıdır. Bu çalışma ‘evren bu mimaride qubit'e kapalıdır’ DEMEZ; eldeki veri böyle bir iddiayı desteklemez. Söylediği şey daha dar ve daha dürüsttür: Γ-noktası ve küçük süper-hücre seviyesinde görülen parazit-moment / qubit-kalitesi ödünleşmesi ölçek değişimi altında kararlı değildir; bu yaklaşım bu ölçekte bu soruyu çözmemektedir.

Bu çalışmanın sınırları açıkça belirtilir:

- Bu ölçekte yalnızca birkaç (~4–6) yapısal-olarak-farklı geometri mevcuttur; ayrık örgüde sürekli bir kiralite/şekil ekseni fiziksel olarak yoktur. Bu, istatistiksel gücü doğal olarak sınırlar.

- Tek kusur (single defect) incelenmiştir; çift-nokta (exchange) rejimi kapsam dışıdır.

- Ortalama-alan Hubbard, momentler için standart fakat yaklaşık bir yöntemdir; Γ-noktasında çözülmüştür (lokalize moment için olağan seçim). Daha yüksek mertebeli korelasyon etkileri ele alınmamıştır.

- Test iki Nd değeriyle (3 ve 4) sınırlıdır. Tek bir ölçek sıçraması verdict'i kesinleştirmeye yeter (sıralama zaten döner), fakat ödünleşmenin ÇOK daha büyük ölçeklerde yeniden ortaya çıkıp çıkmayacağı bu çalışmanın kapsamı dışındadır.

- Nd=4'te ~1700–2000 atom, bu yaklaşımın masaüstü-uygulanabilir sınırına ulaşıldığını gösterir; daha büyük ölçek farklı bir hesap altyapısı gerektirir. Bu, sonucun kendisi kadar önemli bir pratik bulgudur.

**Yine de, dikkat çekici bir gözlem dürüstçe kaydedilmelidir:** Nd=3'te |dN|'in büyümesiyle q ve p'nin BİRLİKTE artması, Vanević'in işaret ettiği ‘kusuru tanımlayan dengesizlik aynı zamanda parazit momenti besler’ hipoteziyle yön olarak tutarlıdır. Bu çalışma bu sinyalin ölçek-dayanıklı OLMADIĞINI gösterir; sinyalin kendisinin sahteliğini değil. Ayrım önemlidir ve gelecek çalışma için bir başlangıç noktasıdır.

## 7. Sonuç

Delikli-grafen kusur-qubitinde parazit kenar momenti ile qubit kalitesi arasındaki ödünleşme sorusu, ilk kez aynı taramada ve önceden-kayıtlı bir dayanıklılık ölçütüyle test edilmiştir. Küçük ölçekte (Nd=3) görülen ödünleşme eğilimi, süper-hücre büyütüldüğünde (Nd=4) işaret olarak tersine döndüğü için, önceden tanımlanmış ölçüt bunu bir ARTIFACT olarak sınıflandırır. Bu negatif sonuç, mevcut literatürün yapmadığı bir şeydir: ‘armchair = bedava’ anlatısını naif kabul etmek yerine, ölçek-dayanıklılığını şeffaf ve önceden-kayıtlı bir ölçütle sınamak ve negatifi raporlamak. Sonucun bilimsel değeri keşifte değil, bir fikrin sonuna kadar dürüstçe izlenmesinde ve aracın verdiği olumsuz cevabın gizlenmeden kabul edilmesindedir.

## Kaynaklar

Pedersen, T. G. ve ark., “Graphene Antidot Lattices: Designed Defects and Spin Qubits”, Phys. Rev. Lett. 100, 136804 (2008). arXiv:0802.4019.

Pedersen grubu, “Large and stable band gaps in spin-polarized graphene antidot lattices”, Phys. Rev. B 88, 195418 (2013).

Vanević, M. ve ark., antidot kenar/midgap yerel momentlerinin dekoherans rolü üzerine (2009). arXiv:0903.0918.

Cirio, M., Brennen, G. K., Twamley, J., levite süperiletken mikroküre qubit önerisi, Phys. Rev. Lett. (2012).

Johnsson, M., Brennen, G. K., Twamley, J., Sci. Rep. (2016).

*Uygulama: tek dosyalık hesaplama hattı (holey_graphene_qubit.py, numpy/scipy/joblib). Tüm sayısal değerler önceden-kayıtlı dayanıklılık testinin doğrudan çıktısıdır; sonradan düzenlenmemiştir.*
