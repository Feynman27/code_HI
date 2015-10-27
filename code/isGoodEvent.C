int isGoodBCID_1(int bcid); // 2010 + 2011
int isGoodBCID_4(int bcid); // 2010
int isGoodBCID_8(int bcid); // 2011
int isGoodBCID_16(int bcid); // 2010
int isGoodBCID_65(int bcid); // 2010
int isGoodBCID_112(int bcid); // 2010
int isGoodBCID_113(int bcid); // 2010
int isGoodBCID_129(int bcid); // 2010
int isGoodBCID_168(int bcid); // 2011
int isGoodBCID_344(int bcid); // 2011
int isGoodBCID_356(int bcid); // 2011
int isGoodBCID_pp(int bcid); // 2010

int isGoodEvent(int run, int lb, int bcid)
{
	switch(run)
	{
// run 2011:

		case 194382:	if(lb>=161 && lb<=253)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194374:	if(lb>=139 && lb<=180)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194370:	if(lb>=227 && lb<=320)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194287:	return 0;

		case 194193:	if(lb>=218 && lb<=587)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194192:	if(lb>=117 && lb<=439)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194179:	if(lb>=102 && lb<=285)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194163:	if(lb>=111 && lb<=518)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194160:	if((lb>=139 && lb<=145) || (lb>=295 && lb<=483))
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194121:	if(lb>=3 && lb<=415)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194079:	return 0;

		case 194061:	if(lb>=334 && lb<=698)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 194060:	if(lb>=3 && lb<=265)
						{
							return isGoodBCID_344(bcid);
						}
						else  return 0;

		case 194059:	return 0;

		case 194056:	return 0;
						
		case 194017:	if(lb>=42 && lb<=376)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193890:	if(lb>=23 && lb<=93)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193834:	if(lb>=21 && lb<=437)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193826:	if(lb>=327 && lb<=694)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193825:	if((lb>=161 && lb<=195) || (lb>=198 && lb<=517))
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193823:	if(lb>=89 && lb<=118)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193795:	if(lb>=673 && lb<=1008)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193718:	if(lb>=126 && lb<=433)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193687:	if(lb>=11 && lb<=18)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193679:	if(lb>=317 && lb<=664)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193662:	if(lb>=109 && lb<=435)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193655:	if(lb>=246 && lb<=526)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193641:	if(lb>=108 && lb<=274)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193604:	if(lb>=166 && lb<=596)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193599:	if(lb>=121 && lb<=137)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;
						
		case 193558:	if(lb>=99 && lb<=364)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193546:	if(lb>=263 && lb<=650)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193494:	if(lb>=3 && lb<=67)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193493:	if((lb>=3 && lb<=95) || (lb>=102 && lb<=104))
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;
						
		case 193492:	if(lb>=5 && lb<=10)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193491:	if(lb>=146 && lb<=303)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193481:	if(lb>=186 && lb<=568)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193463:	if(lb>=78 && lb<=386)
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;
						
		case 193447:	if((lb>=131 && lb<=150) || (lb>=190 && lb<=219) || (lb>=224 && lb<=479))
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193412:	if((lb>=3 && lb<=10) || (lb>=43 && lb<=386))
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193411:	return 0;
						
		case 193403:	if((lb>=5 && lb<=402) || (lb>=404 && lb<=457))
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193402:	return 0;

		case 193321:	if((lb>=39 && lb<=72) || (lb>=87 && lb<=366))
						{
							return isGoodBCID_356(bcid);
						}
						else  return 0;

		case 193295:	if((lb>=4 && lb<=7) || (lb>=36 && lb<=58))
						{
							return isGoodBCID_168(bcid);
						}
						else  return 0;

		case 193294:	return 0;

		case 193291:	if((lb>=442 && lb<=660) || (lb>=704 && lb<=705))
						{
							return isGoodBCID_168(bcid);
						}
						else  return 0;
						
		case 193270:	if(lb>=47 && lb<=460)
						{
							return isGoodBCID_8(bcid);
						}
						else  return 0;
						
		case 193211:	if(lb>=672 && lb<=931)
						{
							return isGoodBCID_1(bcid);
						}
						else  return 0;	
																											
		default: return -1;

	
// run 2010:

		case 168595:	return 0;
						
		case 168603:	return 0;
						
		case 168665:	if((lb>=256 && lb<=439) || (lb>=501 && lb<=559) || (lb>=561 && lb<=580))
						{
							return isGoodBCID_1(bcid);
						}
						else  return 0;
						
		case 168726:	if(lb>=4 && lb<=30)
						{
							return isGoodBCID_1(bcid);
						}
						else  return 0;
						
		case 168759:	if(lb>=227 && lb<=505)
						{
							return isGoodBCID_4(bcid);
						}
						else  return 0;
						
		case 168875:	return 0;
						
		case 169045:	if((lb>=70 && lb<=189) || (lb>=192 && lb<=404))
						{
							return isGoodBCID_65(bcid);
						}
						else  return 0;
						
		case 169136:	if(lb>=248 && lb<=792)
						{
							return isGoodBCID_65(bcid);
						}
						else  return 0;
						
		case 169175:	if(lb>=298 && lb<=755)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;
						
		case 169206:	if(lb>=327 && lb<=700)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;
						
		case 169207:	if(lb>=400 && lb<=455)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;				
						
		case 169223:	if((lb>=89 && lb<=187) || (lb>=191 && lb<=409) || (lb>=440 && lb<=510))
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;
		
		case 169224:	if(lb>=11 && lb<=105)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;	
						
		case 169226:	if((lb>=171 && lb<=177) || (lb>=179 && lb<=508) || (lb>=514 && lb<=540))
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;		
							
		case 169270:	if(lb>=124 && lb<=434)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;				
			
		case 169564:	if((lb>=396 && lb<=500) || (lb>=511 && lb<=535))
						{
							return isGoodBCID_112(bcid);
						}
						else  return 0;		
									
		case 169566:	if(lb>=7 && lb<=198)
						{
							return isGoodBCID_112(bcid);
						}
						else  return 0;
						
		case 169567:	if(lb>=227 && lb<=340)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;
						
		case 169627:	if((lb>=272 && lb<=306) || (lb>=308 && lb<=657))
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 169648:	if(lb>=71 && lb<=133)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;
						
		case 169693:	if((lb>=149 && lb<=152) || (lb>=160 && lb<=723))
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;
		
		case 169750:	if(lb>=281 && lb<=648)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 169751:	if(lb>=287 && lb<=449)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;															

		case 169765:	if((lb>=166 && lb<=267) || (lb>=269 && lb<=473))
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 169783:	if(lb>=109 && lb<=171)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 169837:	return 0;

		case 169838:	return 0;						

		case 169839:	if(lb>=3 && lb<=224)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 169864:	if(lb>=220 && lb<=459)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;		
										
		case 169866:	return 0;

		case 169884:	return 0;

		case 169927:	if((lb>=254 && lb<=312) || (lb>=314 && lb<=715) || (lb>=717 && lb<=830))
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;
						
		case 169961:	if(lb>=83 && lb<=408)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 169964:	if(lb>=129 && lb<=220)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 169966:	if(lb>=616 && lb<=736)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 170002:	if((lb>=144 && lb<=244) || (lb>=246 && lb<=269) || (lb>=271 && lb<=568) || (lb>=570 && lb<=608))
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 170004:	if((lb>=303 && lb<=325) || (lb>=328 && lb<=359) || (lb>=402 && lb<=525) || (lb>=566 && lb<=795))
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 170015:	if((lb>=415  && lb<=545) || (lb>=556 && lb<=584))
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 170016:	if(lb>=3 && lb<=208)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 170080:	if(lb>=215 && lb<=240)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;	
												
		case 170082:	if(lb>=4 && lb<=120)
						{
							return isGoodBCID_113(bcid);
						}
						else  return 0;

		case 170398:	if((lb>=271 && lb<=341) || (lb>=344 && lb<=652))
						{
							return isGoodBCID_129(bcid);
						}
						else  return 0;

		case 170432:	if(lb>=247 && lb<=252)
						{
							return isGoodBCID_129(bcid);
						}
						else  return 0;

		case 170459:	if((lb>=212 && lb<=250) || (lb>=266 && lb<=506))
						{
							return isGoodBCID_129(bcid);
						}
						else  return 0;	
						
		case 170467:	if(lb>=264 && lb<=584)
						{
							return isGoodBCID_129(bcid);
						}
						else  return 0;	

		case 170482:	if(lb>=139 && lb<=349)
						{
							return isGoodBCID_129(bcid);
						}
						else  return 0;

	//2010 pp@2.76TeV
	case 178264:	if((lb>=85 && lb<=183)||(lb>=401 && lb<=423)||(lb>=425 && lb<=556))
						{
							return isGoodBCID_pp(bcid);
						}
						else  return 0;

	case 178229:	if((lb>=164 && lb<=220)||(lb>=222 && lb<=595)||(lb>=597 && lb<=1062)||(lb>=1092 && lb<=1159)||(lb>=1244 && lb<=1440))

						{
							return isGoodBCID_pp(bcid);
						}
						else  return 0;

	case 178211:	if(lb>=160 && lb<=307)						
						{
							return isGoodBCID_pp(bcid);
						}
						else  return 0;

	case 178163:	if((lb>=680 && lb<=694)||(lb>=709 && lb<=799))

						{
							return isGoodBCID_pp(bcid);
						}
						else  return 0;

	}

}

int isGoodBCID_1(int bcid)
{
	switch(bcid)
	{
		case 1:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_4(int bcid)
{
	switch(bcid)
	{
		case 1:	
		case 892:
		case 1783:
		case 2674:	
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_8(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 401:
		case 892:
		case 1292:
		case 1783:
		case 2183:
		case 2674:
		case 3074:
					return 1;
		default:	return 0;
	}
}
int isGoodBCID_16(int bcid)
{
	switch(bcid)
	{
		case 201:
		case 221:
		case 241:
		case 261:
		case 1092:
		case 1112:
		case 1132:
		case 1152:
		case 1983:
		case 2003:
		case 2023:
		case 2043:
		case 2874:
		case 2894:
		case 2914:
		case 2934:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_65(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 107:
		case 127:
		case 147:
		case 167:
		case 273:
		case 293:
		case 313:
		case 333:
		case 439:
		case 459:
		case 479:
		case 499:
		case 832:
		case 852:
		case 872:
		case 892:
		case 998:
		case 1018:
		case 1038:
		case 1058:
		case 1164:
		case 1184:
		case 1204:
		case 1224:
		case 1330:
		case 1350:
		case 1370:
		case 1390:
		case 1496:
		case 1516:
		case 1536:
		case 1556:
		case 1723:
		case 1743:
		case 1763:
		case 1783:
		case 1889:
		case 1909:
		case 1929:
		case 1949:
		case 2055:
		case 2075:
		case 2095:
		case 2115:
		case 2221:
		case 2241:
		case 2261:
		case 2281:
		case 2614:
		case 2634:
		case 2654:
		case 2674:
		case 2780:
		case 2800:
		case 2820:
		case 2840:
		case 2946:
		case 2966:
		case 2986:
		case 3006:
		case 3112:
		case 3132:
		case 3152:
		case 3172:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_112(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 40:
		case 60:
		case 80:
		case 100:
		case 139:
		case 159:
		case 179:
		case 199:
		case 238:
		case 258:
		case 278:
		case 298:
		case 337:
		case 357:
		case 377:
		case 397:
		case 436:
		case 456:
		case 476:
		case 496:
		case 733:
		case 753:
		case 773:
		case 793:
		case 832:
		case 852:
		case 872:
		case 892:
		case 931:
		case 951:
		case 971:
		case 991:
		case 1030:
		case 1050:
		case 1070:
		case 1090:
		case 1129:
		case 1149:
		case 1169:
		case 1189:
		case 1228:
		case 1248:
		case 1268:
		case 1288:
		case 1327:
		case 1367:
		case 1387:
		case 1426:
		case 1446:
		case 1466:
		case 1486:
		case 1525:
		case 1545:
		case 1565:
		case 1585:
		case 1624:
		case 1644:
		case 1664:
		case 1684:
		case 1723:
		case 1743:
		case 1763:
		case 1783:
		case 1822:
		case 1842:
		case 1862:
		case 1882:
		case 1921:
		case 1941:
		case 1961:
		case 1981:
		case 2020:
		case 2040:
		case 2060:
		case 2080:
		case 2119:
		case 2139:
		case 2159:
		case 2179:
		case 2218:
		case 2238:
		case 2258:
		case 2278:
		case 2515:
		case 2535:
		case 2555:
		case 2575:
		case 2614:
		case 2634:
		case 2654:
		case 2674:
		case 2713:
		case 2733:
		case 2753:
		case 2773:
		case 2812:
		case 2832:
		case 2852:
		case 2872:
		case 2911:
		case 2931:
		case 2951:
		case 2971:
		case 3010:
		case 3030:
		case 3050:
		case 3070:
		case 3109:
		case 3129:
		case 3149:
		case 3169:	
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_113(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 40:
		case 60:
		case 80:
		case 100:
		case 139:
		case 159:
		case 179:
		case 199:
		case 238:
		case 258:
		case 278:
		case 298:
		case 337:
		case 357:
		case 377:
		case 397:
		case 436:
		case 456:
		case 476:
		case 496:
		case 733:
		case 753:
		case 773:
		case 793:
		case 832:
		case 852:
		case 872:
		case 892:
		case 931:
		case 951:
		case 971:
		case 991:
		case 1030:
		case 1050:
		case 1070:
		case 1090:
		case 1129:
		case 1149:
		case 1169:
		case 1189:
		case 1228:
		case 1248:
		case 1268:
		case 1288:
		case 1327:
		case 1347:
		case 1367:
		case 1387:
		case 1426:
		case 1446:
		case 1466:
		case 1486:
		case 1525:
		case 1545:
		case 1565:
		case 1585:
		case 1624:
		case 1644:
		case 1664:
		case 1684:
		case 1723:
		case 1743:
		case 1763:
		case 1783:
		case 1822:
		case 1842:
		case 1862:
		case 1882:
		case 1921:
		case 1941:
		case 1961:
		case 1981:
		case 2020:
		case 2040:
		case 2060:
		case 2080:
		case 2119:
		case 2139:
		case 2159:
		case 2179:
		case 2218:
		case 2238:
		case 2258:
		case 2278:
		case 2515:
		case 2535:
		case 2555:
		case 2575:
		case 2614:
		case 2634:
		case 2654:
		case 2674:
		case 2713:
		case 2733:
		case 2753:
		case 2773:
		case 2812:
		case 2832:
		case 2852:
		case 2872:
		case 2911:
		case 2931:
		case 2951:
		case 2971:
		case 3010:
		case 3030:
		case 3050:
		case 3070:
		case 3109:
		case 3129:
		case 3149:
		case 3169:	
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_129(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 40:
		case 60:
		case 80:
		case 100:
		case 120:
		case 140:
		case 160:
		case 180:
		case 218:
		case 238:
		case 258:
		case 278:
		case 298:
		case 318:
		case 338:
		case 358:
		case 396:
		case 416:
		case 436:
		case 456:
		case 476:
		case 496:
		case 516:
		case 536:
		case 752:
		case 772:
		case 792:
		case 812:
		case 832:
		case 852:
		case 872:
		case 892:
		case 931:
		case 951:
		case 971:
		case 991:
		case 1011:
		case 1031:
		case 1051:
		case 1071:
		case 1109:
		case 1129:
		case 1149:
		case 1169:
		case 1189:
		case 1209:
		case 1229:
		case 1249:
		case 1287:
		case 1307:
		case 1327:
		case 1347:
		case 1367:
		case 1387:
		case 1407:
		case 1427:
		case 1465:
		case 1485:
		case 1505:
		case 1525:
		case 1545:
		case 1565:
		case 1585:
		case 1605:
		case 1643:
		case 1663:
		case 1683:
		case 1703:
		case 1723:
		case 1743:
		case 1763:
		case 1783:
		case 1822:
		case 1842:
		case 1862:
		case 1882:
		case 1902:
		case 1922:
		case 1942:
		case 1962:
		case 2000:
		case 2020:
		case 2040:
		case 2060:
		case 2080:
		case 2100:
		case 2120:
		case 2140:
		case 2178:
		case 2198:
		case 2218:
		case 2238:
		case 2258:
		case 2278:
		case 2298:
		case 2318:
		case 2534:
		case 2554:
		case 2574:
		case 2594:
		case 2614:
		case 2634:
		case 2654:
		case 2674:
		case 2713:
		case 2733:
		case 2753:
		case 2773:
		case 2793:
		case 2813:
		case 2833:
		case 2853:
		case 2891:
		case 2911:
		case 2931:
		case 2951:
		case 2971:
		case 2991:
		case 3011:
		case 3031:
		case 3069:
		case 3089:
		case 3109:
		case 3129:
		case 3149:
		case 3169:
		case 3189:
		case 3209:	
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_168(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 9:
		case 17:
		case 25:
		case 33:
		case 41:
		case 49:
		case 57:
		case 65:
		case 73:
		case 81:
		case 89:
		case 97:
		case 105:
		case 113:
		case 121:
		case 129:
		case 137:
		case 145:
		case 153:
		case 161:
		case 169:
		case 177:
		case 185:
		case 443:
		case 451:
		case 459:
		case 467:
		case 475:
		case 483:
		case 491:
		case 499:
		case 507:
		case 515:
		case 523:
		case 531:
		case 539:
		case 547:
		case 555:
		case 563:
		case 571:
		case 579:
		case 587:
		case 595:
		case 603:
		case 611:
		case 619:
		case 627:
		case 892:
		case 900:
		case 908:
		case 916:
		case 924:
		case 932:
		case 940:
		case 948:
		case 956:
		case 964:
		case 972:
		case 980:
		case 988:
		case 996:
		case 1004:
		case 1012:
		case 1020:
		case 1028:
		case 1036:
		case 1044:
		case 1052:
		case 1060:
		case 1068:
		case 1076:
		case 1334:
		case 1342:
		case 1350:
		case 1358:
		case 1366:
		case 1374:
		case 1382:
		case 1390:
		case 1398:
		case 1406:
		case 1414:
		case 1422:
		case 1430:
		case 1438:
		case 1446:
		case 1454:
		case 1462:
		case 1470:
		case 1478:
		case 1486:
		case 1494:
		case 1502:
		case 1510:
		case 1518:
		case 1783:
		case 1791:
		case 1799:
		case 1807:
		case 1815:
		case 1823:
		case 1831:
		case 1839:
		case 1847:
		case 1855:
		case 1863:
		case 1871:
		case 1879:
		case 1887:
		case 1895:
		case 1903:
		case 1911:
		case 1919:
		case 1927:
		case 1935:
		case 1943:
		case 1951:
		case 1959:
		case 1967:
		case 2225:
		case 2233:
		case 2241:
		case 2249:
		case 2257:
		case 2265:
		case 2273:
		case 2281:
		case 2289:
		case 2297:
		case 2305:
		case 2313:
		case 2321:
		case 2329:
		case 2337:
		case 2345:
		case 2353:
		case 2361:
		case 2369:
		case 2377:
		case 2385:
		case 2393:
		case 2401:
		case 2409:
		case 2674:
		case 2682:
		case 2690:
		case 2698:
		case 2706:
		case 2714:
		case 2722:
		case 2730:
		case 2738:
		case 2746:
		case 2754:
		case 2762:
		case 2770:
		case 2778:
		case 2786:
		case 2794:
		case 2802:
		case 2810:
		case 2818:
		case 2826:
		case 2834:
		case 2842:
		case 2850:
		case 2858:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_344(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 9:
		case 17:
		case 25:
		case 33:
		case 41:
		case 49:
		case 57:
		case 65:
		case 73:
		case 81:
		case 89:
		case 97:
		case 105:
		case 113:
		case 121:
		case 129:
		case 137:
		case 145:
		case 153:
		case 161:
		case 169:
		case 177:
		case 185:
		case 222:
		case 230:
		case 238:
		case 246:
		case 254:
		case 262:
		case 270:
		case 278:
		case 286:
		case 294:
		case 302:
		case 310:
		case 318:
		case 326:
		case 334:
		case 342:
		case 350:
		case 358:
		case 366:
		case 374:
		case 382:
		case 390:
		case 398:
		case 406:
		case 443:
		case 451:
		case 459:
		case 467:
		case 475:
		case 483:
		case 491:
		case 499:
		case 507:
		case 515:
		case 523:
		case 531:
		case 539:
		case 547:
		case 555:
		case 563:
		case 571:
		case 579:
		case 680:
		case 688:
		case 696:
		case 704:
		case 712:
		case 720:
		case 728:
		case 736:
		case 744:
		case 752:
		case 760:
		case 768:
		case 776:
		case 784:
		case 792:
		case 800:
		case 808:
		case 816:
		case 824:
		case 832:
		case 840:
		case 848:
		case 892:
		case 900:
		case 908:
		case 916:
		case 924:
		case 932:
		case 940:
		case 948:
		case 956:
		case 964:
		case 972:
		case 980:
		case 988:
		case 996:
		case 1004:
		case 1012:
		case 1020:
		case 1028:
		case 1036:
		case 1044:
		case 1052:
		case 1060:
		case 1068:
		case 1076:
		case 1113:
		case 1121:
		case 1129:
		case 1137:
		case 1145:
		case 1153:
		case 1161:
		case 1169:
		case 1177:
		case 1185:
		case 1193:
		case 1201:
		case 1209:
		case 1217:
		case 1225:
		case 1233:
		case 1241:
		case 1249:
		case 1257:
		case 1265:
		case 1273:
		case 1281:
		case 1289:
		case 1297:
		case 1334:
		case 1342:
		case 1350:
		case 1358:
		case 1366:
		case 1374:
		case 1382:
		case 1390:
		case 1398:
		case 1406:
		case 1414:
		case 1422:
		case 1430:
		case 1438:
		case 1446:
		case 1454:
		case 1462:
		case 1470:
		case 1478:
		case 1486:
		case 1494:
		case 1502:
		case 1510:
		case 1518:
		case 1555:
		case 1563:
		case 1571:
		case 1579:
		case 1587:
		case 1595:
		case 1603:
		case 1611:
		case 1619:
		case 1627:
		case 1635:
		case 1643:
		case 1651:
		case 1659:
		case 1667:
		case 1675:
		case 1683:
		case 1691:
		case 1699:
		case 1707:
		case 1715:
		case 1723:
		case 1731:
		case 1739:
		case 1783:
		case 1791:
		case 1799:
		case 1807:
		case 1815:
		case 1823:
		case 1831:
		case 1839:
		case 1847:
		case 1855:
		case 1863:
		case 1871:
		case 1879:
		case 1887:
		case 1895:
		case 1903:
		case 1911:
		case 1919:
		case 1927:
		case 1935:
		case 1943:
		case 1951:
		case 1959:
		case 1967:
		case 2004:
		case 2012:
		case 2020:
		case 2028:
		case 2036:
		case 2044:
		case 2052:
		case 2060:
		case 2068:
		case 2076:
		case 2084:
		case 2092:
		case 2100:
		case 2108:
		case 2116:
		case 2124:
		case 2132:
		case 2140:
		case 2148:
		case 2156:
		case 2164:
		case 2172:
		case 2180:
		case 2188:
		case 2225:
		case 2233:
		case 2241:
		case 2249:
		case 2257:
		case 2265:
		case 2273:
		case 2281:
		case 2289:
		case 2297:
		case 2305:
		case 2313:
		case 2321:
		case 2329:
		case 2337:
		case 2345:
		case 2353:
		case 2361:
		case 2369:
		case 2377:
		case 2385:
		case 2393:
		case 2401:
		case 2409:
		case 2462:
		case 2470:
		case 2478:
		case 2486:
		case 2494:
		case 2502:
		case 2510:
		case 2518:
		case 2526:
		case 2534:
		case 2542:
		case 2550:
		case 2558:
		case 2566:
		case 2574:
		case 2582:
		case 2590:
		case 2598:
		case 2606:
		case 2614:
		case 2622:
		case 2630:
		case 2674:
		case 2682:
		case 2690:
		case 2698:
		case 2706:
		case 2714:
		case 2722:
		case 2730:
		case 2738:
		case 2746:
		case 2754:
		case 2762:
		case 2770:
		case 2778:
		case 2786:
		case 2794:
		case 2802:
		case 2810:
		case 2818:
		case 2826:
		case 2834:
		case 2842:
		case 2850:
		case 2858:
		case 2895:
		case 2903:
		case 2911:
		case 2919:
		case 2927:
		case 2935:
		case 2943:
		case 2951:
		case 2959:
		case 2967:
		case 2975:
		case 2983:
		case 2991:
		case 2999:
		case 3007:
		case 3015:
		case 3023:
		case 3031:
		case 3039:
		case 3047:
		case 3055:
		case 3063:
		case 3071:
		case 3079:
		case 3116:
		case 3124:
		case 3132:
		case 3140:
		case 3148:
		case 3156:
		case 3164:
		case 3172:
		case 3180:
		case 3188:
		case 3196:
		case 3204:
		case 3212:
		case 3220:
		case 3228:
		case 3236:
		case 3244:
		case 3252:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_pp(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 38:
		case 75:
		case 112:
		case 412:
		case 433:
		case 454:
		case 475:
		case 496:
		case 517:
		case 538:
		case 559:
		case 722:
		case 743:
		case 764:
		case 785:
		case 806:
		case 827:
		case 848:
		case 869:
		case 1081:
		case 1102:
		case 1123:
		case 1144:
		case 1165:
		case 1186:
		case 1207:
		case 1228:
		case 1303:
		case 1324:
		case 1345:
		case 1366:
		case 1387:
		case 1408:
		case 1429:
		case 1450:
		case 2194:
		case 2215:
		case 2236:
		case 2257:
		case 2278:
		case 2299:
		case 2320:
		case 2341:
		case 2507:
		case 2528:
		case 2549:
		case 2570:
		case 2591:
		case 2612:
		case 2633:
		case 2654:
		case 2866:
		case 2887:
		case 2908:
		case 2929:
		case 2950:
		case 2971:
		case 2992:
		case 3013:
		case 3085:
		case 3106:
		case 3127:
		case 3148:
		case 3169:
		case 3190:
		case 3211:
		case 3232:
		return 1;
		default: return 0;
	}

}

int isGoodBCID_356(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 9:
		case 17:
		case 25:
		case 33:
		case 41:
		case 49:
		case 57:
		case 65:
		case 73:
		case 81:
		case 89:
		case 97:
		case 105:
		case 113:
		case 121:
		case 129:
		case 137:
		case 145:
		case 153:
		case 161:
		case 169:
		case 177:
		case 185:
		case 222:
		case 230:
		case 238:
		case 246:
		case 254:
		case 262:
		case 270:
		case 278:
		case 286:
		case 294:
		case 302:
		case 310:
		case 318:
		case 326:
		case 334:
		case 342:
		case 350:
		case 358:
		case 366:
		case 374:
		case 382:
		case 390:
		case 398:
		case 406:
		case 443:
		case 451:
		case 459:
		case 467:
		case 475:
		case 483:
		case 491:
		case 499:
		case 507:
		case 515:
		case 523:
		case 531:
		case 539:
		case 547:
		case 555:
		case 563:
		case 571:
		case 579:
		case 587:
		case 595:
		case 603:
		case 611:
		case 619:
		case 627:
		case 680:
		case 688:
		case 696:
		case 704:
		case 712:
		case 720:
		case 728:
		case 736:
		case 744:
		case 752:
		case 760:
		case 768:
		case 776:
		case 784:
		case 792:
		case 800:
		case 808:
		case 816:
		case 824:
		case 832:
		case 840:
		case 848:
		case 892:
		case 900:
		case 908:
		case 916:
		case 924:
		case 932:
		case 940:
		case 948:
		case 956:
		case 964:
		case 972:
		case 980:
		case 988:
		case 996:
		case 1004:
		case 1012:
		case 1020:
		case 1028:
		case 1036:
		case 1044:
		case 1052:
		case 1060:
		case 1068:
		case 1076:
		case 1113:
		case 1121:
		case 1129:
		case 1137:
		case 1145:
		case 1153:
		case 1161:
		case 1169:
		case 1177:
		case 1185:
		case 1193:
		case 1201:
		case 1209:
		case 1217:
		case 1225:
		case 1233:
		case 1241:
		case 1249:
		case 1257:
		case 1265:
		case 1273:
		case 1281:
		case 1289:
		case 1297:
		case 1334:
		case 1342:
		case 1350:
		case 1358:
		case 1366:
		case 1374:
		case 1382:
		case 1390:
		case 1398:
		case 1406:
		case 1414:
		case 1422:
		case 1430:
		case 1438:
		case 1446:
		case 1454:
		case 1462:
		case 1470:
		case 1478:
		case 1486:
		case 1494:
		case 1502:
		case 1510:
		case 1518:
		case 1555:
		case 1563:
		case 1571:
		case 1579:
		case 1587:
		case 1595:
		case 1603:
		case 1611:
		case 1619:
		case 1627:
		case 1635:
		case 1643:
		case 1651:
		case 1659:
		case 1667:
		case 1675:
		case 1683:
		case 1691:
		case 1699:
		case 1707:
		case 1715:
		case 1723:
		case 1731:
		case 1739:
		case 1783:
		case 1791:
		case 1799:
		case 1807:
		case 1815:
		case 1823:
		case 1831:
		case 1839:
		case 1847:
		case 1855:
		case 1863:
		case 1871:
		case 1879:
		case 1887:
		case 1895:
		case 1903:
		case 1911:
		case 1919:
		case 1927:
		case 1935:
		case 1943:
		case 1951:
		case 1959:
		case 1967:
		case 2004:
		case 2012:
		case 2020:
		case 2028:
		case 2036:
		case 2044:
		case 2052:
		case 2060:
		case 2068:
		case 2076:
		case 2084:
		case 2092:
		case 2100:
		case 2108:
		case 2116:
		case 2124:
		case 2132:
		case 2140:
		case 2148:
		case 2156:
		case 2164:
		case 2172:
		case 2180:
		case 2188:
		case 2225:
		case 2233:
		case 2241:
		case 2249:
		case 2257:
		case 2265:
		case 2273:
		case 2281:
		case 2289:
		case 2297:
		case 2305:
		case 2313:
		case 2321:
		case 2329:
		case 2337:
		case 2345:
		case 2353:
		case 2361:
		case 2369:
		case 2377:
		case 2385:
		case 2393:
		case 2401:
		case 2409:
		case 2462:
		case 2470:
		case 2478:
		case 2486:
		case 2494:
		case 2502:
		case 2510:
		case 2518:
		case 2526:
		case 2534:
		case 2542:
		case 2550:
		case 2558:
		case 2566:
		case 2574:
		case 2582:
		case 2590:
		case 2598:
		case 2606:
		case 2614:
		case 2622:
		case 2630:
		case 2674:
		case 2682:
		case 2690:
		case 2698:
		case 2706:
		case 2714:
		case 2722:
		case 2730:
		case 2738:
		case 2746:
		case 2754:
		case 2762:
		case 2770:
		case 2778:
		case 2786:
		case 2794:
		case 2802:
		case 2810:
		case 2818:
		case 2826:
		case 2834:
		case 2842:
		case 2850:
		case 2858:
		case 2895:
		case 2903:
		case 2911:
		case 2919:
		case 2927:
		case 2935:
		case 2943:
		case 2951:
		case 2959:
		case 2967:
		case 2975:
		case 2983:
		case 2991:
		case 2999:
		case 3007:
		case 3015:
		case 3023:
		case 3031:
		case 3039:
		case 3047:
		case 3055:
		case 3063:
		case 3071:
		case 3079:
		case 3116:
		case 3124:
		case 3132:
		case 3140:
		case 3148:
		case 3156:
		case 3164:
		case 3172:
		case 3180:
		case 3188:
		case 3196:
		case 3204:
		case 3212:
		case 3220:
		case 3228:
		case 3236:
		case 3244:
		case 3252:
		case 3260:
		case 3268:
		case 3276:
		case 3284:
		case 3292:
		case 3300:
					return 1;
		default:	return 0;
	}
}





