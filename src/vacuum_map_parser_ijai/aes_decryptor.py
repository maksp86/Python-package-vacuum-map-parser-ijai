from Crypto.Cipher import AES
from Crypto.Hash import MD5
from Crypto.Util.Padding import pad, unpad
import base64
import logging

_LOGGER = logging.getLogger(__name__)

isEncryptKeyTypeHex = True

def aes_encrypt(data, key: str):
    cipher = AES.new(key.encode("utf-8"), AES.MODE_ECB)

    encryptedData = cipher.encrypt(
        pad(data.encode("utf-8"), AES.block_size, 'pkcs7'))
    encryptedBase64Str = base64.b64encode(encryptedData).decode("utf-8")

    return encryptedBase64Str

def aes_decrypt(data, key: str):
    parsedKey = key.encode("utf-8")
    if isEncryptKeyTypeHex:
        parsedKey = bytes.fromhex(key)

    cipher = AES.new(parsedKey, AES.MODE_ECB)

    decryptedBytes = cipher.decrypt(data)

    decryptedData = unpad(decryptedBytes, AES.block_size, 'pkcs7')
    
    return bytes.fromhex(decryptedData.decode("utf-8"))


def md5key(string: str, model: str, device_mac: str):
    pjstr = "".join(device_mac.lower().split(":"))

    tempModel = model.split('.')[-1]

    if len(tempModel) == 2:
        tempModel = "00" + tempModel
    elif len(tempModel) == 3:
        tempModel = "0" + tempModel
    elif len(tempModel) > 4:
        tempModel = tempModel[-4:]
        
    tempKey = pjstr + tempModel
    aeskey = aes_encrypt(string, tempKey)

    temp = MD5.new(aeskey.encode('utf-8')).hexdigest()
    if isEncryptKeyTypeHex:
        return temp
    else:
        return temp[8:-8].upper()


def gen_md5_key(wifi_info_sn: str, owner_id: str, device_id: str, model: str, device_mac: str):
    arr = [wifi_info_sn, owner_id, device_id]
    tempString = '+'.join(arr)
    return md5key(tempString, model, device_mac)


def decrypt(data: bytes, wifi_info_sn: str, owner_id: str, device_id: str, model: str, device_mac: str) -> bytes:
    try:
        data = base64.b64decode(data)
    except:
        pass
    return aes_decrypt(data, gen_md5_key(wifi_info_sn, owner_id, device_id, model, device_mac))

